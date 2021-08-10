
if (!require(seqinr, quietly = T)) install.packages("seqinr")
library(seqinr)


args = commandArgs(trailingOnly=TRUE)
if (length(args) >= 1 & (args[1] == TRUE | args[1] == "TRUE")){
	print("Sampling new NZ sequences")
	targetCounties = c("New Zealand")
	n.target = 100

	sampling.df = data.frame(size = c("large"),
						method = c("time"),
						n.world = c(420),
						n.china = 80,
						col = "red")
									

}else{
	print("Sampling sequences for 8 countries ")
	targetCounties = c("Colombia")
	n.target = 250


	sampling.df = data.frame(size = c("small", "large", "small", "large"),
						method = c("time", "time", "active", "active"),
						n.world = c(160, 420),
						n.china = c(40, 80),
						col = c(	rgb(0.8, 0.8, 0, 0.5), 
									rgb(0, 0.8, 0.8, 0.5), 
									rgb(0.8, 0, 0.8, 0.5), 
									rgb(0.8, 0.8, 0.8, 0.5)))

}

set.seed(12345)






AUX_SAMPLE = FALSE


# n.aus = c(14, 34),
num.p = 0
date.p = 1
num_bad_seq = 0

#setwd("E:/STOREDFILES/Research/Marsden2019/CoV/COVID19/sequences/subsampling")



get.metadata = function(filename, seqfile, NZ = FALSE) {


	metadata = readLines(filename)[-1]
	headers = strsplit(readLines(filename)[1], "\t")[[1]]
	
	metadata = metadata[grep(".+", metadata)]


	# Meta data
	#sapply(strsplit(metadata, "\t"), function(ele) length(ele))
	metadata.df = data.frame(strain = sapply(strsplit(metadata, "\t"), function(ele) ele[headers == "strain"]),
							date =  sapply(strsplit(metadata, "\t"), function(ele) ele[headers == "date"]),
							region = sapply(strsplit(metadata, "\t"), function(ele) ele[headers == "region"]),
							country = sapply(strsplit(metadata, "\t"), function(ele) ele[headers == "country"]),
							division = sapply(strsplit(metadata, "\t"), function(ele) ele[headers == "division"]),
							host = sapply(strsplit(metadata, "\t"), function(ele) ele[headers == "host"]),
							sex = sapply(strsplit(metadata, "\t"), function(ele) ele[headers == "sex"]),
							seq = "")
							

	metadata.df$strain = as.character(metadata.df$strain)
	metadata.df$date = as.character(metadata.df$date)
	metadata.df$region = as.character(metadata.df$region)
	metadata.df$country = as.character(metadata.df$country)
	metadata.df$division = as.character(metadata.df$division)
	metadata.df$sex = as.character(metadata.df$sex)
	metadata.df$host = as.character(metadata.df$host)
	metadata.df$seq = as.character(metadata.df$seq)
	if (!NZ) {
		metadata.df$clade = as.character(sapply(strsplit(metadata, "\t"), function(ele) ele[headers == "pangolin_lineage"]))
	}else{
		metadata.df$clade = "?"
	}


	# Remove non latin characters
	metadata.df$strain = iconv(metadata.df$strain, "latin1", "ASCII", sub="")
	metadata.df$country = iconv(metadata.df$country, "latin1", "ASCII", sub="")
	metadata.df$region = iconv(metadata.df$region, "latin1", "ASCII", sub="")
	metadata.df$division = iconv(metadata.df$division, "latin1", "ASCII", sub="")

	
	

	# Correct for inconsistent gender labelling
	metadata.df$sex = 	ifelse(metadata.df$sex == "Male", "M",  
						ifelse(metadata.df$sex == "Female", "F", "?"))
						

	# Just human hosts
	metadata.df = metadata.df[metadata.df$host == "Human",]



	# Only if date has day resolution (YYYY-MM-DD)
	metadata.df = metadata.df[grep("20[0-9][0-9]-[0-9][0-9]-[0-9][0-9]", metadata.df$date),]


	# England, Northern Ireland, Scotland, Wales -> United Kingdom (according to gisaid)
	metadata.df$country = ifelse(metadata.df$country == "England" | 
								metadata.df$country == "Scotland" | 
								metadata.df$country == "Wales" |
								metadata.df$country == "NorthernIreland", "United Kingdom", metadata.df$country)
								
								
	# Hong kong -> china (according to gisaid)
	#metadata.df$country = ifelse(metadata.df$country == "Hong Kong", "China",  metadata.df$country)


	# Remove NZ from world pool
	if (!NZ) {
		metadata.df = metadata.df[metadata.df$country != "New Zealand",]
	}
	


	# Date check
	badDates = which(as.Date(metadata.df$date) > Sys.Date())
	if (length(badDates) > 0){
		print(paste("Removing", metadata.df[badDates,"strain"], "because the date is",  metadata.df[badDates,"date"]))
		metadata.df = metadata.df[-badDates,]
	}
	




	# Remove outliers
	outliers = readLines("https://raw.githubusercontent.com/nextstrain/ncov/master/defaults/exclude.txt")
	for (o in outliers) {
		if (o == "" | substring(o, 1, 1) == "#") {
			next
		}
		if (sum(metadata.df$strain == o) > 0) {
			print(paste("Removing", o))
		}
		metadata.df = metadata.df[metadata.df$strain != o,]
	}



	# Sequences
	fasta = readLines(seqfile)
	acc.lines = grep(">", fasta)
	accessions = gsub(">", "", fasta[acc.lines])
	sequences = character(length(acc.lines))
	for (k in 1:length(acc.lines)) {
		start = acc.lines[k]+1
		end = ifelse(k == length(acc.lines), length(fasta), acc.lines[k+1]-1)
		#print(paste(start, end))
		seq = toupper(gsub(" ", "", paste(fasta[start:end], collapse = "")))
		sequences[k] = seq
	}
	print("sequences loaded")


	# Put sequences in data frame
	matches = 0
	strains = metadata.df$strain
	if (NZ){
		strains = sapply(strsplit(strains, "/"), function(ele) ele[2])
	}
	for (i in 1:length(accessions)) {
		acc = accessions[i]
		
		if (NZ) {
			acc2  = strsplit(acc, "_")[[1]][1]
		}else{
			acc2 = gsub("hCoV-19/", "", acc)
			acc2 = strsplit(acc2, "[|]")[[1]][1]
		}
	

		match = which(strains == acc2)
		if (length(match) > 1) match = match[1]
		
		if (length(match) == 1) {
			
			#seq = toupper(paste(fasta[[acc]], collapse = ""))
			seq = sequences[i]
			metadata.df[match,"seq"] = seq
			matches = matches + 1
		}else{
			if (NZ) print(paste("Cannot find", acc))
		}
		
	}		

	print(paste(matches, "matches"))


	# Remove any sequences without matches
	metadata.df = metadata.df[metadata.df$seq != "",]





	# Correct name of NZ sequences
	metadata.df[metadata.df$country == "New Zealand","strain"] = as.character(ifelse(sapply(metadata.df[metadata.df$country == "New Zealand","strain"], function(ele) length(strsplit(ele, "/")[[1]]) > 1),
						metadata.df[metadata.df$country == "New Zealand","strain"],
						paste("NewZealand", metadata.df[metadata.df$country == "New Zealand","strain"], metadata.df[metadata.df$country == "New Zealand","date"], sep = "/")))




	# Remove sequences with more than 10% ambiguous characters
	badseq = (nchar(metadata.df$seq) - nchar(gsub("N", "", metadata.df$seq))) / nchar(metadata.df$seq) > 0.0025
	metadata.df = metadata.df[!badseq | metadata.df$country == "New Zealand" | metadata.df$country == "Taiwan",]
	num_bad_seq <<- num_bad_seq + sum(badseq)


	metadata.df

}

setwd("~/Dropbox/COVID19/COVID-SouthAmerica/Sequences/rawsequences1/ESR")
metadata.esr.df = get.metadata("metadata.tsv", "sequences.fasta", TRUE)

### Taking the Global sampling available in GISAID 
setwd("~/Dropbox/COVID19/COVID-SouthAmerica/Sequences/rawsequences1/GISAID_Updated/ncov_global")
metadata.global = get.metadata("ncov_global.tsv", "ncov_global.fasta")

### Taking the Global sampling available in GISAID 
#setwd("~/Dropbox/COVID19/COVID-SouthAmerica/Sequences/rawsequences1/Exploratory analysis")
#metadata.global = get.metadata("metadata_04092020.tsv", "sequences_04092020.fasta")

### Taking the total of South American sequences
setwd("~/Dropbox/COVID19/COVID-SouthAmerica/Sequences/rawsequences1/GISAID_Updated")
metadata.SA = get.metadata("SouthAmerica24022021.tsv", "SouthAmericanSequences24022021.fasta")
metadata.df = rbind(metadata.global, metadata.esr.df, metadata.SA)





# Infections per country

setwd("~/Dropbox/COVID19/COVID-SouthAmerica/Sequences/rawsequences1/GISAID_Updated")
country.file = readLines("gisaid_31_march.html")
country.file = country.file[grep("<h5>", country.file)]
countries.infected = gsub("[*]", "", gsub("<.+", "", gsub('.+d6d6d6">', "", country.file)))
num.infected = as.numeric(gsub(",", "", gsub("</strong>.+", "", gsub(".+<strong>", "", country.file))))
infections.df = data.frame(country = countries.infected, num = num.infected, weight = 0)
infections.df$country = as.character(infections.df$country)


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
	
	
	# There is no DRC in the webscrape
	metadata.df = metadata.df[metadata.df$country != "Democratic Republic of the Congo",]
	metadata.df = metadata.df[metadata.df$country != "Belize",]
	metadata.df = metadata.df[metadata.df$country != "Mali",]
	
		

	# Any missing countries?
	cnt.remove = character(0)
	for (i in 1:nrow(metadata.df)){
		cnt1 = metadata.df[i,"country"]
		match = which(df$country == cnt1)
		if (length(match) == 0){
			cnt.remove = c(cnt.remove, cnt1)
			print(paste("Cannot find", cnt1, "in df"))
		}else{
			#print(paste(cnt1, "=", df$country[match][1]))
		}
	}
	for (cnt1 in cnt.remove){
		metadata.df = metadata.df[metadata.df$country != cnt1,]
	}
	
	
	earliestDate = getDateFromColname(colnames(df)[2])
	latestDate = getDateFromColname(colnames(df)[ncol(df)])

	#df = df[order(df$country),]
	worldometer[[x]] = df

}



numCols = min(ncol(worldometer$death), ncol(worldometer$recov), ncol(worldometer$infect))
worldometer$death = worldometer$death[,1:numCols]
worldometer$recov = worldometer$recov[,1:numCols]
worldometer$infect = worldometer$infect[,1:numCols]
if (!all(colnames(worldometer$death) == colnames(worldometer$recov)) ||
	!all(colnames(worldometer$death) == colnames(worldometer$infect))){
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
	

	# Cumulative -> non-cumulative
	#df[,3:ncol(df)] = df[,3:ncol(df)] - df[,2:(ncol(df)-1)]
	

	worldometer$active = df

}


print(paste("Pulling from bumbeishvili/covid19-daily-data/. Latest case date:", latestDate))

##




sample.rows.by.method = function(method, df, match, nsamples, inft = FALSE){

	if (method == "time") res = sample.rows.by.date(df, match, nsamples, inft)
	if (method == "infected") res = sample.rows.by.infections(df, match, nsamples, inft)
	if (method == "active") res = sample.rows.by.infections.time(df, match, nsamples)
	if (method == "latest") res = sample.rows.by.latest(df, match, nsamples)
	res
}



sample.rows.by.latest = function(df, match, nsamples){

	df$sampled = FALSE
	df$match = FALSE
	df[match,"match"] = TRUE

	dates = as.Date(df[match,"date"])
	cutoff = sort(dates, decreasing = TRUE)[nsamples]
	meetCutoff = df$date >= cutoff
	samples = numeric(0)
	for (i in 1:nsamples){
		
		latest = which(df$sampled == FALSE & df$match == TRUE & meetCutoff)
		latest = latest[length(latest)]
		df[latest,"sampled"] = TRUE
		samples = c(samples, latest)
	
	}
	samples
	
	

}


sample.rows.by.infections.time = function(df, match, nsamples){


	df$sampled = FALSE
	df$match = FALSE
	df[match,"match"] = TRUE
	df$country.date = as.character(apply(df[,c("country", "date")], 1, function(row) paste(row, collapse = "*")))
	infection.dates = as.Date(sapply(colnames(worldometer$active)[-1], function(ele) 
				as.character(getDateFromColname(ele))))
	
	
	# Prob vector
	country.dates = sort(unique(df[df$sampled == FALSE & df$match == TRUE,"country.date"]))


	getProbVector = function(country.dates, quietly = TRUE) {
	
	
		as.numeric(sapply(country.dates, function(cnt.date) {

				cnt = strsplit(cnt.date, "[*]")[[1]][1]
				d8 = as.Date(strsplit(cnt.date, "[*]")[[1]][2])
				dateMatch = which(infection.dates == d8)
				if (length(dateMatch) == 0){
					closestDate = infection.dates[abs(infection.dates - d8) == min(abs(infection.dates - d8))][1]
					if (!quietly) print(paste("Cannot find any infections on", d8, "in", cnt, "will use", closestDate))
					dateMatch = 1
				}
				
				
				
				#x = worldometer$active[worldometer$active$country == cnt, dateMatch + 1]
				
				# Number of new cases on this date
				if (dateMatch == 1) {
					x = worldometer$infect[worldometer$infect$country == cnt, 2]
				}else{
					x = worldometer$infect[worldometer$infect$country == cnt, dateMatch + 2] - worldometer$infect[worldometer$infect$country == cnt, dateMatch + 1]
				}
				#print(x)
				x = sum(x)
				x = max(x, 1)
				x

			}
		))
	
	}
	probs = getProbVector(country.dates, TRUE)

	samples = numeric(0)
	while (length(samples) < nsamples) {
	
		
	
		# Sample a countrydate
		countrydate = sample(country.dates, 1, prob = probs)
		#print(countrydate)
		cnt = strsplit(countrydate, "[*]")[[1]][1]
		day = strsplit(countrydate, "[*]")[[1]][2]
		
		# Randomly select a sequence from that country date
		candidates = which(df$date == day & df$country == cnt & df$sampled == FALSE & df$match == TRUE)
		
		if (length(candidates) != 0) {
		
			rowNum = ifelse(length(candidates) == 1, candidates, sample(candidates, 1))
			samples = c(samples, rowNum)
			
			
			# Flag it as sampled
			df[rowNum,"sampled"] = TRUE
			

			#print(paste(length(samples), countrydate ))
		}
	}
	
	samples

}



sample.rows.by.date = function(df, match, nsamples, inft = FALSE){

	df$sampled = FALSE
	df$match = FALSE
	df[match,"match"] = TRUE
	days = sort(unique(df[df$sampled == FALSE & df$match == TRUE,"date"]))
	samples = numeric(0)
	for (i in 1:nsamples){
	
		# Sample a day
		day = sample(days, 1)
		
		
		# Sample a sequence from that day
		candidates = which(df$date == day & df$sampled == FALSE & df$match == TRUE)
		probs = rep(1, length(candidates))
		
		# Sample each sequence proportional to the number of infections from that country
		if (inft) {
			
			for (j in 1:length(candidates)){
				#cnt = df[candidates[j],"country"]
				numInfected = df[candidates[j],"numInfected"]
				#numCnt = sum(df[df$match == TRUE & df$country == cnt,])
				probs[j] = numInfected# / numCnt
			}
			
		}
		rowNum = ifelse(length(candidates) == 1, candidates, sample(candidates, 1, replace = F, prob = probs))
		samples = c(samples, rowNum)
		
		# Flag it as sampled
		df[rowNum,"sampled"] = TRUE
		
		# Update the days
		days = sort(unique(df[df$sampled == FALSE & df$match == TRUE,"date"]))

	}
	samples
	

}


sample.rows.by.infections = function(df, match, nsamples, inft = FALSE){

	df$sampled = FALSE
	df$match = FALSE
	df[match,"match"] = TRUE
	cnts = unique(df[df$sampled == FALSE & df$match == TRUE,"country"])
	numInfec = as.numeric(sapply(cnts, function(ele) df[df$sampled == FALSE & df$match == TRUE & df$country == ele,"numInfected"][1]))
	samples = numeric(0)
	numInfec[is.na(numInfec)] = round(mean(numInfec[!is.na(numInfec)]))
	for (i in 1:nsamples){
	
		# Sample a country proportional to infection num
		cnt = sample(cnts, 1, replace = F, prob = numInfec)
		
		
		# Sample a sequence from that country
		candidates = which(df$country == cnt & df$sampled == FALSE & df$match == TRUE)
		rowNum = ifelse(length(candidates) == 1, candidates, sample(candidates, 1))
		samples = c(samples, rowNum)
		
		# Flag it as sampled
		df[rowNum,"sampled"] = TRUE
		
		# Update the countries
		cnts = unique(df[df$sampled == FALSE & df$match == TRUE,"country"])
		numInfec = as.numeric(sapply(cnts, function(ele) df[df$sampled == FALSE & df$match == TRUE & df$country == ele,"numInfected"][1]))
		numInfec[is.na(numInfec)] = round(mean(numInfec[!is.na(numInfec)]))

	}
	samples
	

}







# Get numinfected
for (i in 1:nrow(metadata.df)){


	match = which(infections.df$country == metadata.df[i,"country"])
	
	if (length(match) != 1) {
		print(paste("No match for", metadata.df[i,"country"]))
		next
	}
	numInfected = infections.df[match,"num"]
	metadata.df[i,"numInfected"] = numInfected

}



# Plot
png("sequencesVsCases.png", width = 1000, height = 1000, res = 200)
numCases = sum(infections.df$num)
numSeqs = nrow(metadata.df)
plot(0, 0, type = "n", xlim = c(0, 0.3), ylim = c(0, 0.3), xlab = "Proportion of global cases", ylab = "Proportion of global 'good' sequences")
for (cnt in unique(metadata.df$country)) {


	y = sum(metadata.df$country == cnt) / numSeqs
	x = metadata.df[metadata.df$country == cnt,"numInfected"][1] / numCases
	
	if (!is.na(x) && length(x) != 0) {

		points(x, y)
		if (x > 0.03 | y > 0.03) {
			text(x + 0.005, y + 0.005, cnt, adj = 0)
		}
		
		#infections.df[infections.df$country == cnt,"weight"] = metadata.df[metadata.df$country == cnt,"weight"][1]
	}
}
dev.off()





print("------REGULAR SAMPLING------")

pdf("sequencesOverTime.pdf", width = 6, height = 6)
plot.dates = function(dates, add, col = "#DB5079", day1 = "2019-12-10") {

	dates = sort(as.Date(dates))
	day1 = as.Date(day1)
	days = dates - day1
	freqs = sapply(days, function(d) sum(days == d))
	
	if (!add) {
		plot(days, freqs, type = "n", xlab  = paste("Days since", day1), 
				ylab = "Number of sequences", xaxs = "i", yaxs = "i",  main = "Num sequences over time", 
				xlim = c(0, max(days)))
	}

	polygon(c(0, days, max(days) + 5), c(0, freqs, 0), col = col)


}



plot.dates(metadata.df$date, F, "white")





for (targetCnt in targetCounties) {


	cnt.df = metadata.df


	# Include the reference sequence
	if (targetCnt == "New Zealand") {
		include.in = character() #readLines("include.txt")
	}else{
		include.in = "NewZealand/20VR0174/2020-02-27"
	}
	include = cnt.df[1,]
	include = include[-1,]
	for (o in include.in) {
		if (o == "" | substring(o, 1, 1) == "#") {
			next
		}
		if (sum(cnt.df$strain == o) > 0) {
			include = rbind(include, cnt.df[cnt.df$strain == o,])
			print(paste("Including", o))
		}
		cnt.df = cnt.df[cnt.df$strain != o,]
	}



	# Sampling methods
	SAMPLES = list()
	for (s in 1:nrow(sampling.df)){

		ele.df = sampling.df[s,]
		
		# Ref seq
		if (targetCnt == "New Zealand") {
			nworld = ele.df$n.world
		}else{
			nworld = ele.df$n.world - 1
		}
		
		
		# Sample (or use all) sequences in the target country
		n.target.eff = min(n.target, sum(cnt.df$country == targetCnt))
		print(paste(n.target.eff, "sequences from", targetCnt))
		#target = sample.rows.by.method(ele.df$method, cnt.df, cnt.df$country == targetCnt, n.target.eff)
		target = sample.rows.by.method("latest", cnt.df, cnt.df$country == targetCnt, n.target.eff)

		# Sample seqs from China
		china = sample.rows.by.method(ele.df$method, cnt.df, cnt.df$country == "China", ele.df$n.china)
		
		
		#pakistan = sample.rows.by.method("latest", cnt.df, cnt.df$country == "Pakistan", ele.df$n.pakistan)
		#aus = sample.rows.by.method("latest", cnt.df, cnt.df$country == "Australia", ele.df$n.aus)
		#hk = sample.rows.by.method("latest", cnt.df, cnt.df$country == "Hong Kong", ele.df$n.HK)
		#swiss = sample.rows.by.method("latest", cnt.df, cnt.df$country == "Switzerland", ele.df$n.switz)
		#ireland = sample.rows.by.method("latest", cnt.df, cnt.df$country == "Ireland", ele.df$n.ireland)
		#pakistan, aus, hk, swiss, 
	
		
		# Sample from the rest of the world
		world =  sample.rows.by.method(ele.df$method, cnt.df, 	cnt.df$country != targetCnt &
																cnt.df$country != "China", nworld)
		#sample.rows.by.infections.time(cnt.df, cnt.df$country != targetCnt & cnt.df$country != "China", 1000)

		# Sampled alignment
		samples = cnt.df[c(china, world, target),]
		
			
		# Include all of the mandatory sequences
		samples = rbind(include, samples)
		
		
		# Filebase
		filebase = paste(gsub(" ", "", targetCnt), ele.df$size, ele.df$method, sep = ".")
		
		
		# Write the table
		write.table(sort(table(samples$country)), paste(filebase, "tsv", sep = "."), quote = F, sep = "\t",
		row.names = F)
		
		


		# Print sampled sequences to fasta
		headers = paste0(">", samples$strain, "|", samples$date, "|", samples$region, "|", samples$country, "|", samples$division, "|", samples$sex)
		headers = gsub(" ", "_", headers)
		write(paste(paste(headers, samples$seq, sep = "\n"), collapse = "\n"), paste(filebase, "fasta", sep = "."))


		#plot.dates(samples$date, T, as.character(ele.df$col))
		
		
		SAMPLES[[s]] = samples
		
		
	}

}



  axis(1)
axis(2)
legend("topleft", c("full", paste(sampling.df$size, sampling.df$method, sep = "." )), 
	col = c("white", as.character(sampling.df$col)), pch = 16)

dev.off()




print(paste(num_bad_seq, "sequences have too many ambiguous characters and were discarded"))



# Aux sample 
if (AUX_SAMPLE) {

	print("------AUX SAMPLING------")

	# Get the unique list of strains in the other alignments
	strainsInAlignment = character(0)
	for (targetCnt in targetCounties) {
		for (s in 1:nrow(sampling.df)){

			ele.df = sampling.df[s,]
			filebase = paste(ele.df$size, ele.df$method, sep = ".")
			#filebase = paste0("../alignments/", gsub(" ", "", targetCnt), "/aligned.", filebase)
			filebase = paste0(gsub(" ", "", targetCnt), ".", filebase)
			seqs = read.fasta(paste0(filebase, ".fasta"))
			accs = sapply(strsplit(names(seqs), "[|]"), function(ele) ele[1])
			strainsInAlignment = c(strainsInAlignment, accs)
			
		}
	
	}
	
	# Sort
	strainsInAlignment = sort(strainsInAlignment)
	
	# Include reference sequence to the beginning
	strainsInAlignment = c("NewZealand/20VR0174/2020-02-27", strainsInAlignment)
	
	# Remove duplicates
	strainsInAlignment = unique(strainsInAlignment)
	
	remaining.df = metadata.df[sapply(metadata.df$strain, function(ele) !any(ele == strainsInAlignment)),]
	
	

	
	sampled = sample.rows.by.method("time", remaining.df, TRUE, n.aux)
	samples = remaining.df[sampled,]
	
	# Print sampled sequences to fasta
	headers = paste0(">", samples$strain, "|", samples$date, "|", samples$region, "|", samples$country, "|", samples$division, "|", samples$sex)
	headers = gsub(" ", "_", headers)
	write(paste(paste(headers, samples$seq, sep = "\n"), collapse = "\n"), "ClockPrior.fasta")

	
	

}

##############################

### Geting Subsampling by Country

Argentina<-metadata.df %>% filter(country == "Argentina")
table(Argentina$country)
Brazil<-metadata.df %>% filter(country == "Brazil")
Chile<-metadata.df %>% filter(country == "Chile")
Colombia<-metadata.df %>% filter(country == "Colombia")
Ecuador<-metadata.df %>% filter(country == "Ecuador")
Peru<-metadata.df %>% filter(country == "Peru")
Uruguay<-metadata.df %>% filter(country == "Uruguay")
Venezuela<-metadata.df %>% filter(country == "Venezuela")
Suriname<-South_America %>% filter(country == "Suriname")

## getting a Data frame of South American's countries

NseqCountry <- c(length(Argentina$country), length(Brazil$country), length(Chile$country), length(Colombia$country), length(Ecuador$country), length(Peru$country), length(Uruguay$country), length(Venezuela$country), length(Suriname$country))
names(NseqCountry) <- c("Argentina", "Brazil", "Chile", "Colombia", "Ecuador", "Peru", "Uruguay", "Venezuela", "Suriname")
NseqCountry <- as.data.frame(NseqCountry)
NseqCountry

# downsampling data, keep at most 10 sequences per country per date
#sub-sampled the COVID sequences to keep at most 10 isolates per target per date.

des.down.Arg <- Argentina %>% group_by(date) %>% do(slice_sample(., if (count(.)$n < 10) {size = count(.)$n} else {size = 10}, replace = FALSE))
des.down.Bra <- Brazil %>% group_by(date) %>% do(slice_sample(., if (count(.)$n < 10) {size = count(.)$n} else {size = 10}, replace = FALSE))
des.down.Chi <- Chile %>% group_by(date) %>% do(slice_sample(., if (count(.)$n < 10) {size = count(.)$n} else {size = 10}, replace = FALSE))
des.down.Ecu <- Ecuador %>% group_by(date) %>% do(slice_sample(., if (count(.)$n < 10) {size = count(.)$n} else {size = 10}, replace = FALSE))
des.down.Col <- Colombia %>% group_by(date) %>% do(slice_sample(., if (count(.)$n < 10) {size = count(.)$n} else {size = 10}, replace = FALSE))
des.down.Per <- Peru %>% group_by(date) %>% do(slice_sample(., if (count(.)$n < 10) {size = count(.)$n} else {size = 10}, replace = FALSE))
des.down.Uru <- Uruguay %>% group_by(date) %>% do(slice_sample(., if (count(.)$n < 10) {size = count(.)$n} else {size = 10}, replace = FALSE))
des.down.Ven <- Venezuela %>% group_by(date) %>% do(slice_sample(., if (count(.)$n < 10) {size = count(.)$n} else {size = 10}, replace = FALSE))
des.down.Sur <- Suriname %>% group_by(date) %>% do(slice_sample(., if (count(.)$n < 10) {size = count(.)$n} else {size = 10}, replace = FALSE))

des.down = rbind(des.down.Arg, des.down.Bra, des.down.Chi, des.down.Ecu, des.down.Col, des.down.Per, des.down.Uru, des.down.Ven, des.down.Sur)
table(des.down$country)

NseqCountryDownsampling <- c(length(des.down.Arg$country), length(des.down.Bra$country), length(des.down.Chi$country), length(des.down.Col$country), length(des.down.Ecu$country), length(des.down.Per$country), length(des.down.Uru$country), length(des.down.Ven$country), length(des.down.Sur$country))
names(NseqCountryDownsampling) <- c("Argentina", "Brazil", "Chile", "Colombia", "Ecuador", "Peru", "Uruguay", "Venezuela", "Suriname")
NseqCountryDownsampling <- as.data.frame(NseqCountryDownsampling)
NseqCountryDownsampling

Downsampling <- cbind(NseqCountry, NseqCountryDownsampling)
Downsampling


#Extract the desired sequences from each country with downsampling

DownSampling <- DNASequences[which(names(DNASequences) %in% des.down$strain)]
length(DownSampling)
# Print sampled sequences to fasta
headers = paste0(">", des.down$strain, "|", des.down$date, "|", des.down$region, "|", des.down$country, "|", des.down$division, "|", des.down$sex)
headers = gsub(" ", "_", headers)
write(paste(paste(headers, des.down$seq, sep = "\n"), collapse = "\n"), "SouthAmericanDownSampling.fasta")
#write.FASTA(DownSampling, "SouthAmericanDownSampling.fasta")
