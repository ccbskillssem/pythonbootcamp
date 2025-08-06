# CCB Python Bootcamp ("PyCamp") 🐍
## A biannual introductory workshop to using Python for bioinformatics

Presented by the [Center for Computational Biology](https://ccb.berkeley.edu/). Please contact `ccbadmin(at)berkeley.edu` for all questions.

--------------

## What is this?
This is a repository containing the notebooks, datasets, and visual assets for CCB's Python Bootcamp.

## Table of contents
- [Lecture notebooks](#lecture-notebooks)
- [Zoom](#zoom)
- [recordings](#zoom-recordings)
- [Slack](#slack)
- [Datasets](#datasets)
- [FAQs](#faqs)
- [Authorship](#authorship)
- [License](#license)

## Zoom
Link to the zoom is [here](https://berkeley.zoom.us/j/92582693181). Lectures are recorded and will be posted on the PyCamp slack when they are ready!

## Zoom recordings 

| Day | Zoom recording | password |
| ---- | ------ | ------ |
| Monday | [recording](https://drive.google.com/file/d/12TBVjZgkwdjDIbcWTaJiAyaDWDXe3wfy/view?usp=drive_link)| 
| Tuesday | [recording](https://berkeley.zoom.us/rec/share/rHZda36bx9XTDifI41WISkDJ6LK36xgC-5JYgmmCM16dm10x9J36aWCTgGhO88nv.iXwmFSGY2-mgo9qw) | g$g?7w2j
| Wednesday | [recording]()|
| Thursday | [recording]() | 
| Friday | [recording]() | 



## Slack

Make sure to sign up for our slack [here](https://join.slack.com/t/ccbpycampsummer2025/shared_invite/zt-3aaypujs6-6JETne5BwaIwMEndgi3DHQ). All PyCamp-related communication will be done through slack.

## Lecture notebooks
Each notebook is intended to be a lesson on a particular unit of Python. 

| Day | Content | AM Notebook | PM Notebook | Cheat sheet | Solutions |
| ---- | ------ | ------ | ------ | ------ | ------ |
| **Day 1** (Monday) | Interacting with Google Colaboratory. Introduction to types, lists, tuples, sets, dictionaries. Indexing and slicing iterables. Built-in functions. Boolean logic and simple control flow. | [Monday AM](https://colab.research.google.com/drive/1UJdWZqEKnJRaCu_6GJCqdJfJUBiUlbWw) | [Monday PM](https://colab.research.google.com/drive/1XQhkPiDeG_KhMd-vUzsXDvJr_9brbAzG) | [Day 1 cheatsheet](https://drive.google.com/file/d/1qdoaHMW_ogV4yU7MHtOnyVxIR0CxWfr2/view?usp=sharing) | [AM solutions](https://colab.research.google.com/drive/17G-G6MwxhVhVDkkHM31Z1tmzmWTWgwZV) [PM solutions](https://colab.research.google.com/drive/1WvbY4q7TB0bUGLK9_Gi_61aFdedANisW) |
| **Day 2** (Tuesday) | Review basic data structures and methods. Writing custom functions. Introduction to external packages (`numpy`). | [Tuesday AM](https://colab.research.google.com/drive/11Yz1uhsKTHvtktVEUdS3KpY9-lLT6Z-U) | [Tuesday PM](https://colab.research.google.com/drive/1M_CMRnWCXPz7kXo0j5IVhMowb2-1Rzg1) | [Day 2 cheatsheet](https://drive.google.com/file/d/1cPx2l9xlnq5eD26J3ePTj4C4k4QURYr-/view?usp=sharing) | [AM solutions](https://colab.research.google.com/drive/1GKGzAWUSA_CdILWPl25fmKva7_65WLUJ) [PM solutions](https://colab.research.google.com/drive/1xirPNkhi4NboRRjdajE_qou1icXSx-qH) | 
| **Day 3** (Wednesday) | Review and cover more detail on arrays. Mini-projects on data cleaning, exploration, and basic visualization with the `badhealth` and `hepatocellular` datasets. | [Wednesday AM](https://colab.research.google.com/drive/12jgfoY5ojMGvHmEvxDy0jLU0X-UqN3fb) | [Wednesday PM](https://colab.research.google.com/drive/18UvE1-CV1PasuNU7sJBdLlUrE_2vv4ia) | [Day 3 cheatsheet](https://drive.google.com/file/d/1s_DL4l23ihlWRFca5E0odreIUEcvXe7z/view?usp=sharing) | [AM solutions]() [PM solutions]()| 
| **Day 4** (Thursday) | Introduction to `pandas` operations: indexing, slicing, querying, filtering, merging. Exploration of the `clinvar` dataset.| [Thursday AM](https://colab.research.google.com/drive/1LeXMhaDLX5kZz_qFH5eKaI1rPP6qWXcD) | [Thursday PM](https://colab.research.google.com/drive/1nqp0rhHyQV6eDBmNRUbr12n2WZITVGbW) | [Day 4 cheatsheet](https://drive.google.com/file/d/1Fc9Obxer6ymy2gGVrVLJJ7SmVniFy8Rj/view?usp=sharing) | [AM solutions]() [PM solutions]()
| **Day 5** (Friday) | Introduction to Machine Learning and running Python locally | [Friday AM](Coming soon!) | [Friday PM](Coming soon!) | [Day 5 cheatsheet]() | [AM+PM solutions]()

## Datasets
All datasets used in these materials are derived from publicly available data.

| Dataset | Source |
| ---- | ------ |
| `badhealth` | Originally from `Rdatasets`, downloaded [here](https://vincentarelbundock.github.io/Rdatasets/).
| `hepatocellular` | Originally from `Rdatasets`, downloaded [here](https://vincentarelbundock.github.io/Rdatasets/).
| `clinvar` | A table of variants from the [NCBI ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) database. Curated by [Andrew Sharo](https://www.andrewsharo.com/). |
| `facs`, `metadata_facs`, `metadata_droplet` | Single-cell count matrix and metadata tables from the [Tabula Muris Senis](https://www.nature.com/articles/s41586-020-2496-1) dataset. Downloaded using the [UCSC Cell Browser](https://cells.ucsc.edu/?ds=tabula-muris-senis). |
| WBCD | [Wisconsin Breast Cancer Dataset](https://archive.ics.uci.edu/dataset/17/breast+cancer+wisconsin+diagnostic)

## FAQs
Find answers [here](https://ccb.berkeley.edu/ccb-bioinformatics-bootcamp-january-2022-faq/).

## Authorship
Monday-Thursday notebooks and curriculum were developed by [Stacy Li](stacy.li) in Summer 2022, with exercises contributed by PyCamp staff in the years following. Pre-2022 notebooks were developed by previous generations of PyCamp staff. Day 5 notebook was developed by Prakruthi Burra (prakruthi_burra@berkeley.edu) and Carmelle Catamura (carmelle@berkeley.edu)

## License
This repository is licensed under the [CC 4.0 license](https://creativecommons.org/licenses/by/4.0/). You are free to share, redistribute, and adapt these notebooks, although we request that you provide a link to this repository.
