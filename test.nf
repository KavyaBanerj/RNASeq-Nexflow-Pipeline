#!/usr/bin/env nextflow

log.info 'Testing Nextflow on Gitpod!'

process testProcess {
    output:
    stdout into test_out

    script:
    """
    echo 'Hello from Nextflow in Gitpod!'
    """
}

test_out.subscribe { println it }
