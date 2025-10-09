#!/usr/bin/awk -f
BEGIN {
	# Configuration - modify these values
	output_file = "input.sp.in"
	directory_prefix = "T_"
	temp_tag = "BETA"
	temp_list = "1.0,0.4,0.25"
	ntemps = split(temp_list, temps, ",")
	prev_temp_tag = "PREVIOUS_TEMP"
	current_temp_tag = "CURRENT_TEMP"
	# Constant tag replacements
	tags["DENS"] = "0.5"
	tags["HUBBARDU"] = "4"
	tags["VEC1"] = "[2,0]"
	tags["VEC2"] = "[0,2]"
	tags["ITERS"] = "3"
}

# Skip BEGIN block for actual processing
FNR == 1 {
	# Read the entire template into memory
	template = ""
}

{
	template = template $0 "\n"
}

END {
	# Process for each value in range
	for (i = 1; i <= ntemps; i++) {
		outdir = directory_prefix temps[i]
		status = system("mkdir " outdir)
		outfile = outdir "/" output_file
		# Start with template
		content = template
		# Replace varying tag
		beta = 1.0 / temps[i]
		gsub(temp_tag, beta, content)
		gsub(current_temp_tag, temps[i], content)
		prev_temp = "zero"
		if (i != 1) {
			prev_temp = "./" directory_prefix temps[i - 1] "/dca_sp.hdf5"
		}
		gsub(prev_temp_tag, prev_temp, content)
		# Replace constant tags
		for (tag in tags) {
			gsub(tag, tags[tag], content)
		}
		# Write output file
		printf("%s", content) > outfile
		close(output_file)
		print "Created " output_file
	}
}
