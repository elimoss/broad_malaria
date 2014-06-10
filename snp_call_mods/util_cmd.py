# util_cmd.py - This gives a main() function that serves as a nice wrapper
#	around other commands and presents the ability to serve up multiple
#	command-line functions from a single python script.
#
# requires python >= 2.5
#
# dpark@broadinstitute.org
# $Id: util_cmd.py 7351 2013-01-22 22:53:06Z dpark $

import os, tempfile, sys, shutil, optparse, logging

log = logging.getLogger()
tmpDir = None

def setup_logger(log_level):
	loglevel = getattr(logging, log_level.upper(), None)
	assert loglevel, "unrecognized log level: %s" % log_level
	log.setLevel(loglevel)
	h = logging.StreamHandler()
	h.setFormatter(logging.Formatter("%(asctime)s - %(module)s:%(lineno)d:%(funcName)s - %(levelname)s - %(message)s"))
	log.addHandler(h)

def script_name():
	return sys.argv[0].split('/')[-1].rsplit('.',1)[0]

def common_opts(parser, optlist=['tmpDir', 'loglevel']):
	for opt in optlist:
		if opt=='loglevel':
			parser.add_option("--loglevel", dest="loglevel", type='choice',
				help="Verboseness of output.  [default: %default]",
				default='DEBUG',
				choices=['DEBUG','INFO','WARNING','ERROR','CRITICAL','EXCEPTION'])
		elif opt=='tmpDir':
			parser.add_option("--tmpDir", dest="tmpDir", type='string',
				help="Directory for temp files.  [default: %default]",
				default=tmpDir)
		elif opt=='tmpDirKeep':
			parser.add_option("--tmpDirKeep",
				action="store_true", dest="tmpDirKeep",
				help="If set, do not delete the tmpDir if an exception occurs while running.",
				default=False)
		else:
			raise Exception("unrecognized option %s" % opt)
	return parser

def main(commands, version, tool_paths, description=''):
	''' commands: a list of 4-tuples containing the following:
			1. name of command (string, no whitespace)
			2. method to call that takes two arguments, args (a list of required
				arguments) and options (an optparse construct), and returns
				the desired exit code
			3. method to call that returns an optparse parser.  we provide
				the name of the command and a version string for convenience.
			4. the number of required arguments for this command.  If None,
				we allow any number.
			If commands contains exactly one member and the name of the
			only command is None, then we get rid of the whole multi-command
			thing and just present the options for that one function.
		version: the version string to provide to the parser methods of each command
		tool_paths: a dict.  we will set the 'tmpDir' value so that your
			commands will have access to a suggested temp directory
		description: a long string to present as a description of your script
			as a whole if the script is run with no arguments
		log_level: the logging log level to set the log instance to
	'''
	tool_paths['tmpDir'] = find_tmpDir()
	tmpDir = tool_paths['tmpDir']
	
	cmdlist = [x[0] for x in commands]
	commands = dict([(x[0],x[1:]) for x in commands])
	
	if len(cmdlist)==1 and cmdlist[0]==None:
		# only one (nameless) command in this script, simplify
		command = None
		parser = commands[command][1]('', version)
		(options, args) = parser.parse_args()
	else:
		# multiple commands available
		if len(sys.argv) <= 1:
			print "Usage: python %s commandname options" % sys.argv[0]
			if description.strip():
				print description
			print "\ncommands:"
			for cmd in cmdlist:
				print "\t%s" % cmd
			print "\nRun a command with no options for help on that command."
			return
	
		command = sys.argv[1]
		assert command in commands, "command '%s' not recognized" % command
		parser = commands[command][1](command, version)
		(options, args) = parser.parse_args(sys.argv[2:])
	
	if len(args)==0:
		parser.print_help()
		return
	assert commands[command][2]==None or len(args)==commands[command][2], "%d required arguments, got %d." % (commands[command][2], len(args))
	
	setup_logger(not hasattr(options, 'loglevel') and 'DEBUG' or options.loglevel)
	log.info("version: " + parser.get_version())
	log.debug("command line parameters (including implicit defaults): %s %s" % (
		' '.join(args),
		' '.join(
			["%s=%s" % (o.get_opt_string(), vars(options)[o.dest])
			for o in parser.option_list if o.dest!=None])))
	
	## this is so that ajb073 can run R and see emma (assuming we're not on local)
	#RLibs = '/home/unix/dpark/.RLibs' -- eh, we don't run this much from calhoun anymore anyway
	#if os.access(RLibs, os.F_OK):
	#	os.environ['R_LIBS'] = RLibs
	
	if hasattr(options, 'tmpDir'):
		''' If this command has a tmpDir option, use that as a base directory
			and create a subdirectory within it which we will then destroy at
			the end of execution.
		'''
		proposed_dir = 'tmp-%s-%s' % (script_name(),command!=None and command or '')
		if 'LSB_JOBID' in os.environ:
			proposed_dir = 'tmp-%s-%s-%s-%s' % (script_name(),command,os.environ['LSB_JOBID'],os.environ['LSB_JOBINDEX'])
		tempfile.tempdir = tempfile.mkdtemp(prefix='%s-'%proposed_dir, dir=options.tmpDir)
		log.debug("using tempDir: %s" % tempfile.tempdir)
		os.environ['TMPDIR'] = tempfile.tempdir		# this is for running R
		try:
			ret = commands[command][0](args, options)
		except:
			if hasattr(options, 'tmpDirKeep') and options.tmpDirKeep and not (tempfile.tempdir.startswith('/tmp') or tempfile.tempdir.startswith('/local')):
				log.exception("Exception occurred while running %s, saving tmpDir at %s" % (command, tempfile.tempdir))
			else:
				shutil.rmtree(tempfile.tempdir)
			raise
		else:
			shutil.rmtree(tempfile.tempdir)
		return ret
	else:
		# otherwise just run the command
		return commands[command][0](args, options)

def find_tmpDir():
	''' This provides a suggested base directory for a temp dir for use in your
		optparse-based tmpDir option.
	'''
	tmpdir = '/tmp'
	if os.access('/local/scratch', os.X_OK | os.W_OK | os.R_OK):
		tmpdir = '/local/scratch'
	if 'LSB_JOBID' in os.environ:
		# this directory often exists for LSF jobs, but not always.
		# for example, if the job is part of a job array, this directory is called
		# something unpredictable and unfindable, so just use /local/scratch
		proposed_dir = '/local/scratch/%s.tmpdir' % os.environ['LSB_JOBID']
		if os.access(proposed_dir, os.X_OK | os.W_OK | os.R_OK):
			tmpdir = proposed_dir
	return tmpdir
