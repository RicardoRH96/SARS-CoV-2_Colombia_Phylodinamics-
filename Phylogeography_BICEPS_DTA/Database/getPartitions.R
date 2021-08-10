
if (!require(seqinr, quietly = T)) install.packages("seqinr")
library(seqinr)
if (!require(ape, quietly = T)) install.packages("ape")
library(ape)
if (!require(rjson , quietly = T)) install.packages("rjson")
library(rjson)


dirichlet_k = 10000

# Usage:
# Rscript getPartitions.R <alignment.fasta> <alignment.nexus>
# where <alignment.fasta>  is the alignment (input)
# and <alignment.nexus> is the file to save alignment and partitions to (output)
args = commandArgs(trailingOnly=TRUE)


CLEAVE_ORF1 = FALSE
CODONS = TRUE

if (FALSE) {
	setwd("~/Dropbox/COVID19/COVIDSouthAmerica/Sequences/rawsequences1/alignment")
	args = c("~/Dropbox/COVID19/COVIDSouthAmerica/Sequences/rawsequences1/alignment/aligned.large.active.fasta")
}


if (length(args) != 1) {
	stop("Use: Rscript getPartitions.R <alignment.fasta>")
}



partitions.df = read.csv("refseqs/partitions.csv", header = T)
fasta_ref = readLines("refseqs/aligned_ref.fasta")
fasta_new_aln = read.fasta(args[1])




for (col in colnames(partitions.df)) partitions.df[,col] = as.character(partitions.df[,col])



# Reference sequence1 (ref alignment)
ref.seq1.gapped = fasta_ref[2]


# Reference sequence2 (ref alignment)
ref.seq2.name = gsub(">", "", gsub("[|].+", "", fasta_ref[3]))
ref.seq2.gapped = fasta_ref[4]


# Reference sequence2 (new alignment)
new.seq2.names = gsub("[|].+", "", names(fasta_new_aln))
new.seq2.line = which(new.seq2.names == ref.seq2.name)
if (length(new.seq2.line) == 0) {
	stop(paste("cannot find reference sequence", ref.seq2.name, "in", args[1]))
}
if (length(new.seq2.line) > 1) {
	stop(paste("found duplicate sequences", ref.seq2.name, "in", args[1]))
}
new.seq2.gapped = paste(fasta_new_aln[names(fasta_new_aln)[new.seq2.line]][[1]], collapse = "")


# Mapping between reference alignment positions (G) and ungapped positions (U) of refseq1 and refseq1
G_U1 = data.frame(G = 1:nchar(ref.seq1.gapped), U = NA)
G_U2 = G_U1
U_pos_1 = 1
U_pos_2 = 1
for (i in 1:nrow(G_U1)){
	G = G_U1[i,"G"]
	
	# Refseq1
	char1 = substring(ref.seq1.gapped, G, G)
	if (char1 != "-") {
		G_U1[G,"U"] = U_pos_1
		U_pos_1 = U_pos_1 + 1
	}
	
	# Refseq2
	char2 = substring(ref.seq2.gapped, G, G)
	if (char2 != "-") {
		G_U2[G,"U"] = U_pos_2
		U_pos_2 = U_pos_2 + 1
	}
	

}

# Mapping between new alignment positions (G) and ungapped positions (U) of the newseq2 (aka refseq2)
G_U3 = data.frame(G = 1:nchar(new.seq2.gapped), U = NA)
U_pos_3 = 1
for (i in 1:nrow(G_U3)){
	G = G_U3[i,"G"]

	# Newseq2
	char3 = substring(new.seq2.gapped, G, G)
	if (char3 != "-") {
		G_U3[G,"U"] = U_pos_3
		U_pos_3 = U_pos_3 + 1
	}
}


#plot(0, 0, type = "n", xlim = c(1,nchar(new.seq2.gapped)), ylim = c(0, 1))

cols = lapply(rep("red", nrow( partitions.df) ), function(ele) ele)
names(cols) = partitions.df$Partition
cols[["5primeUTR"]] = "white"
cols[["3primeUTR"]] = "white"
cols[["ORF1a"]] = "#696969"
cols[["ORF1b"]] = "#d3d3d3"

cols[["Spike"]] = "#ffdf00"
cols[["ORF3a"]] = "#f4adb4"

cols[["Envelope"]] = "orange"
cols[["Membrane"]] = "brown"
cols[["ORF6"]] = "#f4adb4"
cols[["ORF7a"]] = "#f4adb4"
cols[["ORF8"]] = "#f4adb4"
cols[["Nucleocapsid"]] = "#98fb98"
cols[["ORF10"]] = "#f4adb4"

cols[["IGR"]] = "white"

covered_alignment_area = rep(FALSE, nchar(new.seq2.gapped))
for (p in partitions.df$Partition) {


	if (p != "IGR") {
		print(p)
		


		# Ungapped position in 1st reference sequence
		start.stop = as.character(partitions.df[partitions.df$Partition == p,"Ref1UngappedPos"])
		start.U1 = as.numeric(strsplit(start.stop, "-")[[1]][1])
		stop.U1 = as.numeric(strsplit(start.stop, "-")[[1]][2])
		
		
		# Gapped position in ref alignment
		start.G = G_U1[which(G_U1$U == start.U1),"G"]
		stop.G = G_U1[which(G_U1$U == stop.U1),"G"]
		
		
		
		# Ungapped position in 2nd reference sequence
		start.U2 = G_U2[which(G_U2$G == start.G),"U"]
		stop.U2 = G_U2[which(G_U2$G == stop.G),"U"]
		partitions.df[partitions.df$Partition == p,"Ref2UngappedPos"] = paste(start.U2, stop.U2, sep = "-")
		
		
		
		# Gapped position in the new alignment
		start.G.new = G_U3[which(G_U3$U == start.U2),"G"]
		stop.G.new = G_U3[which(G_U3$U == stop.U2),"G"]
		if (p == "5primeRUBBISH") start.G.new = 1
		partitions.df[partitions.df$Partition == p,"AlnGappedPos"] = paste(start.G.new, stop.G.new, sep = "-")
		
		
		# For finding intergenic regions
		covered_alignment_area[start.G.new:stop.G.new] = TRUE
		
		
		# Ignore poly(A) tail
		if (p == "3primeRUBBISH") covered_alignment_area[stop.G.new:nchar(new.seq2.gapped)] = TRUE
		
		
				
		
		#rect(start.G.new, 0.2, stop.G.new+1, 0.8, col = cols[[p]])
		#text((start.G.new + stop.G.new) / 2, 0.5, p)
		
		
		

	}

}





# Intergenic regions
IGRs = paste(which(covered_alignment_area == FALSE), collapse = ",")
partitions.df[partitions.df$Partition == "IGR","AlnGappedPos"] = IGRs


#write.csv(partitions.df, "refseqs/partitions.csv", quote = F, row.names = F)


if (CLEAVE_ORF1) {

	# Expand ORF1a and 1b into its 15 composite proteins
	polyprotein.df = read.csv("refseqs/polyprotein1ab.csv", header = T)
	fasta_ptn_aln = readLines("refseqs/polyprotein1ab_ref.fasta")

	for (col in colnames(polyprotein.df)) polyprotein.df[,col] = as.character(polyprotein.df[,col])

	# Protein reference sequence1 (ref alignment)
	ref.seq1.ptn.gapped = fasta_ptn_aln[2]


	# Protein reference sequence2 (ref alignment)
	ref.seq2.ptn.gapped = fasta_ptn_aln[4]

			
			
			
	# Mapping between ref protein aligment positions (G) and ungapped positions (U) of protein refseq1 and 2
	G_U_ptn1 = data.frame(G = 1:nchar(ref.seq1.ptn.gapped), U = NA)
	G_U_ptn2 = G_U_ptn1
	U_ptn_pos_1 = 1
	U_ptn_pos_2 = 1
	for (i in 1:nrow(G_U_ptn1)){
		G = G_U_ptn1[i,"G"]

		# Protein seq 1
		char1 = substring(ref.seq1.ptn.gapped, G, G)
		if (char1 != "-") {
			G_U_ptn1[G,"U"] = U_ptn_pos_1
			U_ptn_pos_1 = U_ptn_pos_1 + 1
		}
		
		# Protein seq 2
		char2 = substring(ref.seq2.ptn.gapped, G, G)
		if (char2 != "-") {
			G_U_ptn2[G,"U"] = U_ptn_pos_2
			U_ptn_pos_2 = U_ptn_pos_2 + 1
		}
		
	}



	# Mapping from polyprotein1ab (protein) to ORF1a + ORF1b (genome)
	ORF1.start.ref = as.numeric(strsplit(partitions.df[partitions.df$Partition == "ORF1a","Ref2UngappedPos"], "-")[[1]][1])
	frameshift.position.ref = as.numeric(strsplit(partitions.df[partitions.df$Partition == "ORF1b","Ref2UngappedPos"], "-")[[1]][1])
	for (p in polyprotein.df$Partition) {

		print(p)
		
		
		# Ungapped position in 1st reference sequence (ptn)
		start.stop = as.character(polyprotein.df[polyprotein.df$Partition == p,"Ref1UngappedPos"])
		start.U1 = as.numeric(strsplit(start.stop, "-")[[1]][1])
		stop.U1 = as.numeric(strsplit(start.stop, "-")[[1]][2])
			
			
		# Gapped position in ref alignment (ptn)
		start.G = G_U_ptn1[which(G_U_ptn1$U == start.U1),"G"]
		stop.G = G_U_ptn1[which(G_U_ptn1$U == stop.U1),"G"]
		

		# Ungapped position in 2nd reference sequence (ptn)
		start.U2.ptn = G_U_ptn2[which(G_U_ptn2$G == start.G),"U"]
		stop.U2.ptn = G_U_ptn2[which(G_U_ptn2$G == stop.G),"U"]
		polyprotein.df[polyprotein.df$Partition == p,"Ref2UngappedPos"] = paste(start.U2.ptn, stop.U2.ptn, sep = "-")



		# Map from polyprotein position to genome position
		start.U2.gen = start.U2.ptn * 3 + ORF1.start.ref - 3
		stop.U2.gen = stop.U2.ptn * 3 + ORF1.start.ref - 1
		
			
		# Correct for the frameshift
		if (start.U2.gen >= frameshift.position.ref) start.U2.gen = start.U2.gen - 1
		if (stop.U2.gen >= frameshift.position.ref) stop.U2.gen = stop.U2.gen - 1

		
		
		# Gapped position in the new alignment
		start.G.new = G_U3[which(G_U3$U == start.U2.gen),"G"]
		stop.G.new = G_U3[which(G_U3$U == stop.U2.gen),"G"]
		if (p == "5primeRUBBISH") start.G.new = 1
		polyprotein.df[polyprotein.df$Partition == p,"AlnGappedPos"] = paste(start.G.new, stop.G.new, sep = "-")


	}
	
	
	#write.csv(polyprotein.df, "refseqs/polyprotein1ab.csv", quote = F, row.names = F)
	
	
	
	# Merge the two data frames
	merged.df = rbind(	partitions.df[partitions.df$Partition == "5primeRUBBISH",],
						partitions.df[partitions.df$Partition == "5primeUTR",],
						polyprotein.df,
						partitions.df[	partitions.df$Partition != "5primeRUBBISH" &
										partitions.df$Partition != "5primeUTR" &
										partitions.df$Partition != "ORF1a" &
										partitions.df$Partition != "ORF1b",])
	
	
	


}else {
	merged.df = partitions.df

}

if (CODONS) {


	# Codon positions
	sites.df = data.frame(site = 1:nchar(new.seq2.gapped), position = 5)
	codons.df = data.frame(Partition = c("codon1", "codon2", "codon3", "noncoding"), AlnGappedPos = "")
	for (col in colnames(codons.df)) codons.df[,col] = as.character(codons.df[,col])


	for (p in partitions.df$Partition) {

		if (p == "5primeRUBBISH" |  p == "3primeRUBBISH") next
		
		
		if (p == "IGR") {
			sites = as.numeric(strsplit(partitions.df[p == partitions.df$Partition,"AlnGappedPos"], ",")[[1]])
			for (site in sites) {

				pos = 4
				sites.df[sites.df$site == site, "position"] = min(sites.df[sites.df$site == site, "position"], pos)
			
			}
		} else {
		

			start.pos = as.numeric(strsplit(partitions.df[p == partitions.df$Partition,"AlnGappedPos"], "-")[[1]][1])
			stop.pos = as.numeric(strsplit(partitions.df[p == partitions.df$Partition,"AlnGappedPos"], "-")[[1]][2])


			if (p == "ORF1b") {
				start.pos = start.pos - 1
			}


			for (site in start.pos:stop.pos) {
			
			
				noncoding = p == "5primeUTR" | p == "3primeUTR" | p == "IGR"
				pos = ifelse(noncoding, 4, (site-start.pos) %% 3 + 1)
				sites.df[sites.df$site == site, "position"] = min(sites.df[sites.df$site == site, "position"], pos)
			
			}
			aa.len = (stop.pos - start.pos + 1) / 3
			print(paste(p, aa.len))
		
		}

	}

	codons.df[codons.df$Partition == "codon1","AlnGappedPos"] = paste(sites.df[sites.df$position == 1,"site"], collapse = ",")
	codons.df[codons.df$Partition == "codon2","AlnGappedPos"] = paste(sites.df[sites.df$position == 2,"site"], collapse = ",")
	codons.df[codons.df$Partition == "codon3","AlnGappedPos"] = paste(sites.df[sites.df$position == 3,"site"], collapse = ",")
	codons.df[codons.df$Partition == "noncoding","AlnGappedPos"] = paste(sites.df[sites.df$position == 4,"site"], collapse = ",")
	
	
	#write.table(codons.df, "refseqs/codons.csv", sep = ";", quote = F, row.names = F)


}
#write.csv(merged.df, "refseqs/merged.csv", quote = F, row.names = F)


# Save partitions to txt file ready to insert in a NEXUS file

toprint = "BEGIN ASSUMPTIONS;"

data.df = merged.df
if (CODONS) data.df = codons.df

for (p in data.df$Partition) {



	toprint = paste0(toprint, "\n\t", "charset ", p, " = ", 
				data.df[data.df$Partition == p,"AlnGappedPos"], ";")

}
toprint = paste0(toprint, "\nEND;")
#write(toprint, "assumptions.txt")



# Write gene boundary assumptions
toprint_gene = "BEGIN ASSUMPTIONS;"


for (p in merged.df$Partition) {
	toprint_gene = paste0(toprint_gene, "\n\t", "charset ", p, " = ", 
				merged.df[merged.df$Partition == p,"AlnGappedPos"], ";")

}
toprint_gene = paste0(toprint_gene, "\nEND;")
#write(toprint_gene, "genePartitions.txt")





# Convert to nexus
nexus.name = gsub(".fasta", ".nexus", args[1])
print(paste("Saving to", nexus.name))

fasta.in = read.fasta(args[1])
write.nexus.data(fasta.in, file=nexus.name, format="dna")


# Add assumptions to nexus
nexus.in = readLines(nexus.name)

nexus.out = paste(paste(nexus.in, collapse = "\n"), toprint, sep = "\n") 

write(nexus.out, nexus.name)



# Convert to xml
xml_template = readLines("refseqs/aligned_template.xml")



# Sequences
seqs = character(0)
seq.template = '\t<sequence id="SEQ_ID" spec="Sequence" taxon="SEQ_ID" totalcount="4" value="SEQ_NT"/>'
for (i in 1:length(fasta.in)) {

	seq_nt = paste(fasta.in[[i]], collapse = "")
	seq_id = names(fasta.in)[i]
	seq = gsub("SEQ_ID", seq_id, seq.template)
	seq = gsub("SEQ_NT", seq_nt, seq)
	seqs = c(seqs, seq)
	
}
xml_template = gsub("INSERT_ALIGNMENT_HERE", paste(seqs, collapse = "\n") ,xml_template)


# Partitions
xml_template = gsub("INSERT_CODON1_FILTER_HERE", data.df[data.df$Partition == "codon1","AlnGappedPos"], xml_template)
xml_template = gsub("INSERT_CODON2_FILTER_HERE", data.df[data.df$Partition == "codon2","AlnGappedPos"], xml_template)
xml_template = gsub("INSERT_CODON3_FILTER_HERE", data.df[data.df$Partition == "codon3","AlnGappedPos"], xml_template)
xml_template = gsub("INSERT_NONCODING_FILTER_HERE", data.df[data.df$Partition == "noncoding","AlnGappedPos"], xml_template)
xml_template = gsub("INSERT_WEIGHTS_HERE", paste(sapply(strsplit(data.df[,"AlnGappedPos"], ","), length), collapse = " "), xml_template)


xml.name = gsub(".fasta", ".xml", args[1])
write(xml_template, xml.name)




# JSON
json.name = gsub(".fasta", ".json", args[1])
print(paste("Saving to", json.name))
JSON = list()


# Trim sequence names
seqs = character(0)
seq.template = '\t<sequence id="SEQ_ID" spec="Sequence" taxon="SEQ_ID" totalcount="4" value="SEQ_NT"/>'
for (i in 1:length(fasta.in)) {

	seq_nt = paste(fasta.in[[i]], collapse = "")
	seq_id = strsplit(names(fasta.in)[i], "[|]")[[1]][1]
	seq = gsub("SEQ_ID", seq_id, seq.template)
	seq = gsub("SEQ_NT", seq_nt, seq)
	seqs = c(seqs, seq)
	
}


# Sequences
JSON$sequences = gsub("\"", "'", paste(seqs, collapse = "\n"))
JSON$sequences = paste("\n", JSON$sequences, "\n")



# Partitions
JSON[["filter-codon1"]] = data.df[data.df$Partition == "codon1","AlnGappedPos"]
JSON[["filter-codon2"]] = data.df[data.df$Partition == "codon2","AlnGappedPos"]
JSON[["filter-codon3"]] = data.df[data.df$Partition == "codon3","AlnGappedPos"]
JSON[["filter-noncoding"]] = data.df[data.df$Partition == "noncoding","AlnGappedPos"]
JSON[["partitionweights"]] = paste(sapply(strsplit(data.df[,"AlnGappedPos"], ","), length), collapse = " ")


# Traits
strains = sapply(strsplit(names(fasta.in), "[|]"), function(ele) ele[1])
countries = sapply(strsplit(names(fasta.in), "[|]"), function(ele) ele[4])
dates = as.Date(sapply(strsplit(names(fasta.in), "[|]"), function(ele) ele[2]))
JSON$datetrait = paste(paste0(strains, "=", sapply(strsplit(names(fasta.in), "[|]"), function(ele) ele[2])), collapse = ",")
JSON$regiontrait = paste(paste0(strains, "=", sapply(strsplit(names(fasta.in), "[|]"), function(ele) ele[3])), collapse = ",")
JSON$countrytrait = paste(paste0(strains, "=", countries), collapse = ",")
JSON$divisiontrait = paste(paste0(strains, "=", sapply(strsplit(names(fasta.in), "[|]"), function(ele) ele[5])), collapse = ",")
JSON$sextrait = paste(paste0(strains, "=", sapply(strsplit(names(fasta.in), "[|]"), function(ele) ele[6])), collapse = ",")



# Deme frequencies
wb.df = read.csv("refseqs/worldbank.csv", header = T)


mapToDeme = function(areas, demefile){




	# Map country to deme
	deme.map = read.table(demefile, sep = "\t", header = T)
	wildcard = which(deme.map$from == "*")
	mapped = character(0)
	for (area in areas){
	
		match = which(deme.map$from == area)
		
		if (length(match) == 0){
		
			# Wildcard match
			if (length(wildcard) == 1){
				to = deme.map[wildcard,"to"]
			}else{
				stop(paste("Cannot map", area, "to a deme in", demefile))
			}
		
		}else if (length(match) > 1){
			stop(paste(area, "is mapped to more than 1 deme in", demefile))
		}else{
			to = deme.map[match,"to"]
		}
		mapped = c(mapped, as.character(to))
	}
	
	
	# Dates
	date_intervals = as.Date(seq(from = earliestCaseDate, to = mostRecentTip, by = 7))
	
	# Get num infections per week in each deme
	demes = as.character(unique(deme.map$to))
	demes.df = data.frame(matrix(0, ncol = length(demes) + 1, nrow = length(date_intervals)))
	colnames(demes.df) =  c("times", demes)
	
	demes.df$times = as.numeric(mostRecentTip - date_intervals) / 365
	
	for (i in 1:length(date_intervals)) {
		d8 = date_intervals[i]
		#dateRange = as.Date(seq(from = d8, to = d8 + 6, by = 1))
		dateRange = d8 + 3
		
		cols = format(dateRange, format = "%m.%d.%y")
		cols = paste0("X", cols)
		cols = gsub("[.]0", ".", cols)
		cols = gsub("X0", "X", cols)
		cols = cols[sapply(cols, function(ele) any(ele == colnames(worldometer$active))) ]
		
		
		
		for (deme in demes){
			
			deme_countries = as.character(deme.map[deme.map$to == deme,"from"])
			if (all(deme_countries == "*")) {
				deme_countries = areas[sapply(areas, function(ele) !any(ele == deme.map$from))]
			}
			deme_countries = unique(gsub("_", " ", deme_countries))
			vals = worldometer$active[sapply(worldometer$active$country, function(ele) any(ele == deme_countries)),cols]
			
			
			demes.df[i,deme] = max(sum(vals), 0) + 1
			
		
		}
	
	}
	
	
	
	# Human population frequencies
	demeFreqs = list()
	for (d in demes) demeFreqs[[d]] = 0
	popsum = 0
	for (cnt in unique(areas)){
	
		cnt2 = gsub("_", " ", cnt)
		if (cnt2 == "USA") cnt2 = "United States"
		if (cnt2 == "South Korea") cnt2 = "Korea, Rep."
		if (cnt2 == "Russia") cnt2 = "Russian Federation"
		if (cnt2 == "Brunei") cnt2 = "Brunei Darussalam"
		if (cnt2 == "Slovakia") cnt2 = "Slovak Republic"
		
		
		
		if (cnt2 == "Taiwan"){
			pop = 23780000
		}else {
			match = which(wb.df$Country.Name == cnt2)
			if (length(match) == 0){
				stop(paste("Cannot find", cnt2, "in worldbank"))
			}
			pop = as.numeric(wb.df[match,"X2018"])
		}
		#print(paste("The population of", cnt2, "in 2018 was", pop/1000000, "million"))
		
		
		# Map to deme
		match = which(deme.map$from == cnt)
		
		if (length(match) == 0){
			to = as.character(deme.map[wildcard,"to"])
		}else{
			to = as.character(deme.map[match,"to"])
		}
		
		#print(paste("The population of", cnt2, "(", to, ")", "in 2018 was", pop/1000000, "million"))
		demeFreqs[[to]] = demeFreqs[[to]] + pop
		popsum = popsum + pop
	
	}
	
	demeFreqsNorm = demeFreqs
	for (d in demes) demeFreqsNorm[[d]] = dirichlet_k * demeFreqsNorm[[d]] / popsum
	

	toReturn = list(mapped = mapped, epochs = demes.df, demeFreqs = demeFreqsNorm, popSizes = demeFreqs)
	toReturn
	
}



# Infections per country per day
infections.day.df = read.csv("https://raw.githubusercontent.com/bumbeishvili/covid19-daily-data/master/time_series_19-covid-Confirmed.csv", header = T)
recoveries.day.df = read.csv("https://raw.githubusercontent.com/bumbeishvili/covid19-daily-data/master/time_series_19-covid-Recovered.csv", header = T)
deaths.day.df = read.csv("https://raw.githubusercontent.com/bumbeishvili/covid19-daily-data/master/time_series_19-covid-Deaths.csv", header = T)



getDateFromColname = function(colname) {
	d = gsub("X", "", colname)
	d_split = strsplit(d, "[.]")[[1]]
	as.Date(paste(paste0("20", d_split[3]), d_split[1], d_split[2], sep = "-"))
}


worldometer = list(infect = infections.day.df, recov = recoveries.day.df, death = deaths.day.df)
for (x in names(worldometer)){

	df = worldometer[[x]]
	colnames(df)[colnames(df) == "Country.Region"] = "country" 
	colnames(df)[1] = "country" 
	df = df[,-c(1,3,4)]
	df$country = as.character(df$country)
	
	# Match country names
	df[df$country == "UK","country"] = "United Kingdom"
	df[df$country == "UAE","country"] = "United Arab Emirates"
	#df[,3:ncol(df)] = df[,3:ncol(df)] - df[,2:(ncol(df)-1)]
	earliestCaseDate = as.Date(getDateFromColname(colnames(df)[2]))


	#df = df[order(df$country),]
	worldometer[[x]] = df

}


numCols = min(ncol(worldometer$death), ncol(worldometer$recov), ncol(worldometer$infect))
worldometer$death = worldometer$death[,1:numCols]
worldometer$recov = worldometer$recov[,1:numCols]
worldometer$infect = worldometer$infect[,1:numCols]

if (!all(colnames(worldometer$death) == colnames(worldometer$recov)) ||
	!all(colnames(worldometer$death) == colnames(worldometer$infect))){
	print(colnames(worldometer$death))
	print(colnames(worldometer$recov))
	print(colnames(worldometer$infect))
	stop("Mismatching worldometer column names")
}


# Active cases
worldometer$active = worldometer$infect 
for (cnt in worldometer$active$country){


	df = worldometer$active
	infects = as.numeric(worldometer$infect[worldometer$infect$country == cnt,2:ncol(df)])
	recoveries = as.numeric(worldometer$recov[worldometer$recov$country == cnt,2:ncol(df)])
	deaths = as.numeric(worldometer$death[worldometer$death$country == cnt,2:ncol(df)])
	active = infects - recoveries - deaths 
	active = sapply(active, function(a) max(a, 0))
	df[df$country == cnt,2:ncol(df)] = active
	

	worldometer$active = df

}





# Times
mostRecentTip = max(dates)
oldestRecentTip = min(dates)
JSON$mostRecentTip = as.character(mostRecentTip)
JSON$oldestRecentTip = as.character(oldestRecentTip)
JSON$timeSinceAlert3 = as.numeric(mostRecentTip - as.Date("2020-04-28")) / 365
JSON$timeSinceAlert4 = as.numeric(mostRecentTip - as.Date("2020-03-26")) / 365
JSON$timeSinceNZBorderClose = as.numeric(mostRecentTip - as.Date("2020-03-20")) / 365
JSON$timeSinceFirstNZSample = as.numeric(mostRecentTip - as.Date("2020-02-27")) / 365
JSON$timeSinceFirstAvailableCase = as.numeric(mostRecentTip - earliestCaseDate) / 365


# Get the cumulative number of cases in this time period in this list of countries
 getCasesFromCountriesOnDate = function(areas, endDate = as.Date("3000-01-01"), startDate = as.Date("0000-01-01"), active = FALSE ){

	dateCols = as.character(sapply(colnames(worldometer$recov)[-1], function(ele) as.character(getDateFromColname(ele))))
	startDateCol = which(dateCols > startDate)[1] + 1
	endDateCol = which(dateCols <= endDate)
	endDateCol = endDateCol[length(endDateCol)] + 1
	
	
	
	total = 0
	for (cnt in unique(areas)){
	
		cnt2 = gsub("_", " ", cnt)
		if (active) {
  		ncases = as.numeric(worldometer$infect[worldometer$infect$country == cnt2,c(startDateCol:endDateCol)])
			nrec = as.numeric(worldometer$recov[worldometer$recov$country == cnt2,c(startDateCol:endDateCol)])
			ndeath = as.numeric(worldometer$death[worldometer$death$country == cnt2,c(startDateCol:endDateCol)])
			#print(ncases - nrec - ndeath)
			ncases = mean(ncases - nrec - ndeath)
			print(paste("There were a daily average of ", ncases, "active cases in", cnt2, "between", startDate, "and", endDate))
		
		}else{
		
			ncases = as.numeric(worldometer$infect[worldometer$infect$country == cnt2,c(startDateCol,endDateCol)])
			ncases = ncases[2] - ncases[1]
			print(paste("There were", ncases, "cases in", cnt2, "between", startDate, "and", endDate))
		}
		

		total = total + ncases

	
	}
	
	
#	 Report actual mean number across period
	if (active) {
		total
	
#Pseudocount of 1
	}else{
		max(c(0, total)) + 1
}
	


}


getDemeThenTime = function(df) {
	df = df[,-1] # Remove time column
	str = rev(c(t(df[,ncol(df):1])))
	str = log(str)
	str = (str - mean(str)) / sd(str)
	paste(str, collapse = " ")
}


# Demes
demeFolderName = strsplit(args[1], "/")[[1]]
demeFolderName = paste0(paste(demeFolderName[-length(demeFolderName)], collapse = "/"), "/demes/")
demeFiles = list.files(path = demeFolderName, pattern = "demes")

for (demeFile in demeFiles) {


	#numDemes = as.numeric(gsub("[.]tsv", "", gsub(".+demes", "", demeFile)))
	cat_name = gsub("[.]tsv", "", demeFile)
	
	traitXdemes = mapToDeme(countries, paste0(demeFolderName, demeFile))
	JSON[[cat_name]] = paste(paste0(strains, "=", traitXdemes$mapped), collapse = ",")
	
	
	# Human frequencies
	freqs = paste(as.numeric(round(unlist(traitXdemes$demeFreqs), 1)), collapse = " ")
	JSON[[paste0(cat_name, "_demeDirichletFreqs")]] = freqs

	# Population sizes (for DST)
	for (col in colnames(traitXdemes$epoch)){
		JSON[[paste0(cat_name, "_", col)]] = paste(rev(traitXdemes$epochs[,col]), collapse = " ")
	}

	# Normalised population sizes (for Mascot)
	JSON[[paste0(cat_name, "_pop")]] = getDemeThenTime(traitXdemes$epochs)
	
	
	# 2 deme model (World and Target?)
	demes = as.character(sort(unique(traitXdemes$mapped)))
	if (demes[1] == "Target" & demes[2] == "World"){
	
	
		
		
		
		# Origin prior
		originOffset = (mostRecentTip - as.Date("2019-12-20"))  / 365
		originS = 0.4
		originMean = log(0.15) - 0.5*originS^2
		
		JSON$originPrior = paste0('<LogNormal id="originPriorDistribution" name="distr" S="', originS, '" M="', originMean, '" offset="', originOffset,'"/>')
	
	
		# Compute upper limit for sampling proportion at 2 dates
		firstTargetSequence = min(dates[traitXdemes$mapped == "Target"])
		firstWorldSequence = min(dates[traitXdemes$mapped == "World"])
		lastTargetSequence = max(dates[traitXdemes$mapped == "Target"])
		lastWorldSequence = max(dates[traitXdemes$mapped == "World"])
		JSON$firstTargetSequence = as.character(firstTargetSequence)
		JSON$firstWorldSequence = as.character(firstWorldSequence)
		JSON$lastTargetSequence = as.character(lastTargetSequence)
		JSON$lastWorldSequence = as.character(lastWorldSequence)

		JSON$timeSinceFirstTargetSequence = as.numeric(mostRecentTip - firstTargetSequence) / 365
		JSON$timeSinceFirstWorldSequence = as.numeric(mostRecentTip - firstWorldSequence) / 365
		daydiff = as.numeric(firstTargetSequence - firstWorldSequence)
		targetCnt = unique(countries[traitXdemes$mapped == "Target"])
		worldCnt = unique(countries[traitXdemes$mapped == "World"])
		if (daydiff < 5) {
			print(paste("First target sequence is", daydiff, "days after the first world sequence. Please revise model."))
		}

		# 1) from 1st TARGET sequence to most recent tip
		nTarget = sum(traitXdemes$mapped == "Target" & dates <= mostRecentTip & dates > firstTargetSequence)
		nWorld = sum(traitXdemes$mapped == "World" & dates <= mostRecentTip & dates > firstTargetSequence)
		numRemovs = getCasesFromCountriesOnDate(targetCnt, mostRecentTip, firstTargetSequence)
		JSON$targetSPUpperTarget = nTarget / (numRemovs)


		# 2) from 1st WORLD sequence to most recent tip
		nWorld = sum(traitXdemes$mapped == "World" & dates <= mostRecentTip & dates > firstWorldSequence)
		numRemovs = getCasesFromCountriesOnDate(worldCnt, mostRecentTip, firstWorldSequence)
		JSON$worldSPUpperWorld = nWorld / (numRemovs)
		
		
		# Average number of active cases from 1 Jan to 20th March
		JSON$TargetAverageActiveCasesPreBorder = getCasesFromCountriesOnDate(targetCnt, "2020-03-20", "2020-01-01", TRUE)
		JSON$WorldAverageActiveCasesPreBorder = getCasesFromCountriesOnDate(worldCnt, "2020-03-20", "2020-01-01", TRUE)
		JSON$TargetAverageActiveCasesPostBorder = getCasesFromCountriesOnDate(targetCnt, "2020-04-30", "2020-03-21", TRUE)
		JSON$WorldAverageActiveCasesPostBorder = getCasesFromCountriesOnDate(worldCnt, "2020-04-30", "2020-03-21", TRUE)		
		
		
		
		# Migration rate priors
		dt_beforeIntro = as.numeric(firstTargetSequence - 2019-11-17) / 365
		dt_afterIntro = as.numeric(lastTargetSequence - firstTargetSequence) / 365
		kappa = 0.5
		k = 0.05
		A = 0
		D = 0
		S = 0.5
		if (targetCnt == "New_Zealand"){
			A = 0.52
		}else if (targetCnt == "Australia"){
			A = 0.64
		}else if (targetCnt == "Iceland"){
			A = 0.19
		}else if (targetCnt == "Taiwan"){
			A = 0.8
		} 
		
		print(paste("Country", targetCnt, A))


		JSON$m_wd_before = 1 / (dt_beforeIntro * sum(traitXdemes$mapped == "World") * kappa)
		JSON$m_dw_before = traitXdemes$popSizes[["World"]] / traitXdemes$popSizes[["Target"]] * JSON$m_wd_before

		JSON$m_wd_after = (A * sum(traitXdemes$mapped == "Target")) / (dt_afterIntro * sum(traitXdemes$mapped == "World") * kappa)
		JSON$m_dw_after = traitXdemes$popSizes[["World"]] / traitXdemes$popSizes[["Target"]] * JSON$m_wd_after

		JSON$m_wd_late = k * JSON$m_wd_after
		JSON$m_dw_late = k * JSON$m_dw_after
		

		
		getM = function(rate){
			log(rate) - 0.5*S*S
		}
		
		if (A > 0){
		
		
			JSON$worldPopulation = traitXdemes$popSizes[["World"]]
			JSON$targetPopulation = traitXdemes$popSizes[["Target"]]
		
			
			JSON$EntryPreCaseMPrior = 
				paste0('<LogNormal id="EntryPreCaseMPrior" name="distr" S="', S, '" M="', getM(JSON$m_wd_before), '"/>')
			JSON$ExitPreCaseMPrior = 
				paste0('<LogNormal id="ExitPreCaseMPrior" name="distr" S="', S, '" M="', getM(JSON$m_dw_before), '"/>')


			JSON$EntryPreBorderMPrior = 
				paste0('<LogNormal id="EntryPreBorderMPrior" name="distr" S="', S, '" M="', getM(JSON$m_wd_after), '"/>')
			JSON$ExitPreBorderMPrior = 
				paste0('<LogNormal id="ExitPreBorderMPrior" name="distr" S="', S, '" M="', getM(JSON$m_dw_after), '"/>')
				

			JSON$EntryPostBorderMPrior = 
				paste0('<LogNormal id="EntryPostBorderMPrior" name="distr" S="', S, '" M="', getM(JSON$m_wd_late), '"/>')
			JSON$ExitPostBorderMPrior = 
				paste0('<LogNormal id="ExitPostBorderMPrior" name="distr" S="', S, '" M="', getM(JSON$m_dw_late), '"/>')
		
		}
		
		
	
	}
	
	
	


}



# Deme frequencies
wb.df = read.csv("refseqs/worldbank.csv", header = T)


x = as.character(rjson::toJSON(JSON, indent=1))
x = gsub("\\\\t", "\t", x)
x = gsub("\\\\n", "\n", x)
write(x, json.name)