#!/bin/bash

nextflow run main.nf -profile aws -c nextflow.config
rm -f .nextflow.log* 