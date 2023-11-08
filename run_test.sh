#!/usr/bin/env bash

nextflow run . -profile test,docker -resume
