import urllib.request
import os

DIR_NAME = "data"

def downloadFiles():
    # Make a new directory for downloadable content
    if not os.path.exists(DIR_NAME):
        os.mkdir(DIR_NAME)
    os.chdir(DIR_NAME)

    # Download the file from `url` and save it locally under `file_name`:
    urllib.request.urlretrieve("https://storage.googleapis.com/kagglesdsdata/datasets/668974/1178022/sequence.gb.txt?X-Goog-Algorithm=GOOG4-RSA-SHA256&X-Goog-Credential=gcp-kaggle-com%40kaggle-161607.iam.gserviceaccount.com%2F20210526%2Fauto%2Fstorage%2Fgoog4_request&X-Goog-Date=20210526T231851Z&X-Goog-Expires=259199&X-Goog-SignedHeaders=host&X-Goog-Signature=38465784d866fd77ea5b50c8bb7292f60f9ae15bc9685d27ad598dbe8230c628cefd800b88c2a1bcc2b4bcdcd690b1aa43cbbed608a6db34e7bc319b6094afa15332858ba7a319a1ea90e2be9f75247bd593a779ca75660fb84b7a2ba0bbedfd20cb7f6afddc89b3658cf7ca500f9acf9834b12f0b12301cc0ef53f0b71588a62ae7c849f0e01c910e92e89f54c148ccfcb10b0ef3bb2191c551aaa9c21533fc4c5790ba509618e1c84590713431387dcf38e0b223ff5545b0a223d5ace1991a9f4f174f1d49ae4872e80b25276898436975616b9da13a12f33dfd4323c8785ad880b84e4831dcd4a919cad88276b74875b84517a2f45afdcd1ff65794e98eae", "SARS-CoV-2.gbk")
    urllib.request.urlretrieve("https://www.ebi.ac.uk/ena/browser/api/fasta/MN908947.3", "SARS-CoV-2.fasta")

    os.chdir("..")