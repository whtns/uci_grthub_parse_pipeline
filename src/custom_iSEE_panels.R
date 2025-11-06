#!/usr/bin/env python

import sys
import os
import pandas
import json

proj_name=sys.argv[1]

tracks_conf_file="/var/www/html/jbrowse/" + proj_name + "/trackList.json"
# tracks_conf_file=sys.argv[1]
# bigwig_dir=sys.argv[2]
# metadata_file = sys.argv[3]

bigwig_dir = "/var/www/html/jbrowse/" + proj_name + "/samples/"

bigwig_files = []
for (dirpath, dirnames, filenames) in os.walk(bigwig_dir):
	for filename in filenames:
		if filename.endswith(".bw"):
			bigwig_files.append("samples/" + filename)
	# bigwig_files.extend(filenames)

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
			"label" : os.path.basename(i).replace(".bw", ""),
			"key" : os.path.basename(i).replace(".bw", ""),
			"scale" : "log"
			}
		)
		
# "plugins" : ["NeatHTMLFeatures", "NeatCanvasFeatures", "HideTrackLabels"],		
bigwig_dict = {
	"tracks" : tracks_list,
	"trackMetadata": {
      "sources": [
           { "type": "csv", "url":  "6_seq_fetal/6_seq_metadata.csv" }
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
  "formatVersion" : "1"
}


with open(tracks_conf_file, 'w') as json_file:
  json.dump(bigwig_dict, json_file, indent = 2, sort_keys = False)
# print(json.dumps(person_dict, indent = 4, sort_keys=True))


