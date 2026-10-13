# Exact commands used (DUG, /data/curtin_mwaeor/src/MWAEoR-Pipeline, branch iono-sub-selection-2026-09-03)

## container
    cd /data/curtin_mwaeor/singularity
    singularity pull hyperdrive_0.8.1-bmetrics_cuda12.5.1-ubuntu24.04.sif \
        docker://mwatelescope/hyperdrive:0.8.1-bmetrics-cuda12.5.1-ubuntu24.04
    singularity exec hyperdrive_0.8.1-bmetrics_cuda12.5.1-ubuntu24.04.sif hyperdrive --version        # mwa_hyperdrive 0.8.1
    singularity exec hyperdrive_0.8.1-bmetrics_cuda12.5.1-ubuntu24.04.sif hyperdrive peel --help | grep -c iono-sub-max   # 4
    singularity exec hyperdrive_0.8.1-bmetrics_cuda12.5.1-ubuntu24.04.sif hyperdrive bmetrics --help | head -2

## pipeline runs (launched 2026-10-13 ~15:00 AWST, all three concurrently)
    export http_proxy=http://proxy.per.dug.com:3128 https_proxy=http://proxy.per.dug.com:3128
    nextflow run main.nf -profile dug -with-tower -c eor0low-ionosel-none.conf  --obsids_suffix=-eor0low20260903 --nothumbnail --nopolcomp --novideo --nopsmetrics --nodelayspec --noplotvisqa -name ionosel-none-1
    nextflow run main.nf -profile dug -with-tower -c eor0low-ionosel-pred4.conf --obsids_suffix=-eor0low20260903 --nothumbnail --nopolcomp --novideo --nopsmetrics --nodelayspec --noplotvisqa -name ionosel-pred4-1
    nextflow run main.nf -profile dug -with-tower -c eor0low-ionosel-meas4.conf --obsids_suffix=-eor0low20260903 --nothumbnail --nopolcomp --novideo --nopsmetrics --nodelayspec --noplotvisqa -name ionosel-meas4-1
    nextflow run main.nf -profile dug -with-tower -c eor0low-ionosel-none.conf  --obsids_suffix=-eor0low20260903-quiet --nothumbnail --nopolcomp --novideo --nopsmetrics --nodelayspec --noplotvisqa -name ionosel-none-quiet-1
    nextflow run main.nf -profile dug -with-tower -c eor0low-ionosel-none.conf  --obsids_suffix=-eor0low20260903-event --nothumbnail --nopolcomp --novideo --nopsmetrics --nodelayspec --noplotvisqa -name ionosel-none-event-1
    nextflow run main.nf -profile dug -with-tower -c eor0low-ionosel-pred4.conf --obsids_suffix=-eor0low20260903-quiet --nothumbnail --nopolcomp --novideo --nopsmetrics --nodelayspec --noplotvisqa -name ionosel-pred4-quiet-1
    nextflow run main.nf -profile dug -with-tower -c eor0low-ionosel-pred4.conf --obsids_suffix=-eor0low20260903-event --nothumbnail --nopolcomp --novideo --nopsmetrics --nodelayspec --noplotvisqa -name ionosel-pred4-event-1
    nextflow run main.nf -profile dug -with-tower -c eor0low-ionosel-meas4.conf --obsids_suffix=-eor0low20260903-quiet --nothumbnail --nopolcomp --novideo --nopsmetrics --nodelayspec --noplotvisqa -name ionosel-meas4-quiet-1
    nextflow run main.nf -profile dug -with-tower -c eor0low-ionosel-meas4.conf --obsids_suffix=-eor0low20260903-event --nothumbnail --nopolcomp --novideo --nopsmetrics --nodelayspec --noplotvisqa -name ionosel-meas4-event-1
