#!/usr/bin/env python

import sys
import os
import pandas
import json

proj_name = sys.argv[1]
metadata = sys.argv[2]

tracks_conf_file="/var/www/html/jbrowse/" + proj_name + "/trackList.json"
# tracks_conf_file=sys.argv[1]
# samples_dir=sys.argv[2]
# metadata_file = sys.argv[3]

samples_dir = "/var/www/html/jbrowse/" + proj_name + "/samples/"

bigwig_files = []
bam_files = []
for (dirpath, dirnames, filenames) in os.walk(samples_dir):
	for filename in filenames:
		if filename.endswith(".bw"):
			bigwig_files.append("samples/" + filename)
		if filename.endswith(".bam"):
			bam_files.append("samples/" + filename)
	# sample_files.extend(filenames)

tracks_list = [
	      {
         "chunkSize" : "20000",
         "urlTemplate" : "seq/{refseq}",
         "type" : "SequenceTrack",
         "label" : "DNA",
         "key" : "DNA",
         "metadata": {
             "Category" : "Reference sequence"
             }
         },
         {
           "chunkSize" : "20000",
           "urlTemplate" : "Homo_sapiens.GRCh38.87.sorted.gff3.gz",
           "storeClass" : "JBrowse/Store/SeqFeature/GFF3Tabix",
           "type" : "CanvasFeatures",
           "key" : "genes",
           "label" : "genes"
         }

]

for i in bigwig_files:
	tracks_list.append(
		{
			"urlTemplate" : i,
			"style" : {
				"className" : "image"
			},
			"compress" : "0",
			"type" : "JBrowse/View/Track/Wiggle/XYPlot",
			"label" : os.path.basename(i),
			"key" : "bigwig",
			"scale" : "log"
			}
		)
		
for i in bam_files:
	tracks_list.append(
	  {
         "storeClass"  : "JBrowse/Store/SeqFeature/BAM",
         "urlTemplate" : i,
         "label" : os.path.basename(i),
         "key" : "bam",
         "type"        : "JBrowse/View/Track/Alignments2"
      }
		)
		
# "plugins" : ["NeatHTMLFeatures", "NeatCanvasFeatures", "HideTrackLabels"],		
samples_dict = {
	"tracks" : tracks_list,
	"trackMetadata": {
      "sources": [
           { "type": "csv", "url":  metadata }
      ]
  },
   "trackSelector" : {
      "renameFacets" : {
         "label" : "Track",
         "key" : "Dataset"
      },
      "type" : "Faceted",
      "escapeHTMLInData" : "false",
      "displayColumns" : [
        "batch",
        "key",
        "label",
      ]
   },
  "formatVersion" : "1",
  "dataset_id" : proj_name
}


with open(tracks_conf_file, 'w') as json_file:
  json.dump(samples_dict, json_file, indent = 2, sort_keys = False)
# print(json.dumps(person_dict, indent = 4, sort_keys=True))


