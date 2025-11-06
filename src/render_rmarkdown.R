
args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
	eval(parse(text = args[[i]]))
}

rmarkdown::render(input = rmarkdown_file, output_file = report, params = list(clone_file = clone_file, prepped_pagoda = prepped_pagoda, numbat_dir = numbat_dir))
