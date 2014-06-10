set.seed(1)

#input = '/Volumes/seq_plasmodium/moss/projects/bam_coi/msp/sent064.10/sent064.10_msp1.txt'
input = '/Volumes/seq_plasmodium/moss/projects/bam_coi/msp/sim_mix/msp1/results_50bp.txt'
data = read.table(input, fill = NA)
filt_data = data[data[,1] > 60, ]
components = dim(filt_data)[2]

for (n in 2:components){	
	print(paste(c("Component:", n - 1), sep = " "))
	clust = kmeans(filt_data[,n][!is.na(filt_data[,n])], 2)
	#print(clust$centers)
	if (abs(clust$centers[1] - clust$centers[2]) > 0.1 ){
		if (n==2){
			print(min(clust$centers))	
		}
		else{
			print(max(clust$centers))
		}
	}
}