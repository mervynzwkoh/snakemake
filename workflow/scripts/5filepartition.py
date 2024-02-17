# from selenium.webdriver.common.keys import Keys
# from selenium import webdriver
# from selenium.webdriver.chrome.service import Service
# from webdriver_manager.chrome import ChromeDriverManager
# from selenium.webdriver.support.ui import WebDriverWait
# from selenium.webdriver.common.by import By
# from selenium.webdriver.support import expected_conditions as EC
# import time
import os
import argparse
import xlsxwriter
import shutil, csv
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument(
    "source", help="The source directory of the MasterList"
)  # must be in same directory as msa file
args = parser.parse_args()

# directory with fasta and vcf files must be named "fastavcfloc" with subdirectories "fasta" and "vcf"
fastavcfloc = args.source + "/fastavcfloc"


# obtain absolute path of MasterList
def get_masterlist_path():
    for filename in os.listdir(args.source):
        if filename.endswith(".csv"):
            return os.path.join(args.source + "/", filename)


masterlist = get_masterlist_path()


# obtain number of partitions from clustering algorithm
def get_partcsv_path():
    for filename in os.listdir(args.source + "/distancearrays/clusterlabels/Tamura"):
        if filename.endswith(".csv"):
            return os.path.join(
                args.source + "/distancearrays/clusterlabels/Tamura/", filename
            )


partitions_csv = get_partcsv_path()

max_integer = None
with open(partitions_csv, "r") as csv_file:
    csv_reader = csv.reader(csv_file)
    for row in csv_reader:
        try:
            integer_value = int(row[1])
        except ValueError:
            # The second column does not contain an integer.
            continue
        except IndexError:
            continue

        if max_integer is None or integer_value > max_integer:
            max_integer = integer_value

    num_partitions = max_integer + 1

"""
# this is needed for seqtrack (seqtrack needs collection dates). works only for ncbi database
def getFromAcc(sample):  # gets ReleaseDate from accession number (eg. OUxxx)
    driver = webdriver.Chrome(service=Service(ChromeDriverManager().install()))
    link = "https://www.ncbi.nlm.nih.gov/nuccore/"
    print(sample)
    newLink = link + sample
    print(newLink)
    driver.get(newLink)
    # time.sleep(15.0) #needs some lagtime since website takes a while to load
    # ERR = [my_elem.get_attribute("href") for my_elem in WebDriverWait(driver, 50).until(EC.visibility_of_all_elements_located((By.CSS_SELECTOR, "a[href^='https://www.ncbi.nlm.nih.gov/sra/']")))][0]
    try:
        element_present = EC.presence_of_element_located(
            (By.CSS_SELECTOR, "span[id^='feature']")
        )
        WebDriverWait(driver, 200).until(element_present)  # waits for up to 200 seconds
    except:
        driver.get(newLink)
        time.sleep(30)
        print("Time exceeded")
    features = driver.find_element(By.CSS_SELECTOR, "span[id^='feature']").text
    rDate = features.split('collection_date="', 1)[1][:-1]
    ERR = [
        my_elem.get_attribute("href")
        for my_elem in WebDriverWait(driver, 50).until(
            EC.visibility_of_all_elements_located(
                (By.CSS_SELECTOR, "a[href^='https://www.ncbi.nlm.nih.gov/sra/']")
            )
        )
    ][0]
    # bp = features.split('\n',1)[0].split('..',1)[1]
    driver.close()
    print(rDate)
    return [rDate, ERR]


# src = ""
# out = ""
# with open("AY.122/AccNumbers2.txt", "r") as f:
#     with open("AY.122/Reads.txt", "a") as f1:
#         for l in f.readlines():
#             print(l)
#             f1.write("\n")
#             f1.write(getFromAcc(l))

# src1 = "seqtrack/AY4.10.samples.txt"
# src2 = "seqtrack/AY9.2.samples.txt"
# src3 = "seqtrack/BA1.7.samples.txt"
"""


# make accessionList excel file for each partition
def makeExcel(acc_path, part):
    samples = []
    with open(partitions_csv, "r") as csv_file:
        csv_reader = csv.reader(csv_file)
        for row in csv_reader:
            if row == []:
                continue
            if int(row[1]) == int(part):
                samples.append(str(row[0]))

    with open(masterlist, "r", errors="ignore") as f:
        wbk = acc_path + "/" + str(part) + "accessionList.xlsx"
        workbook = xlsxwriter.Workbook(wbk)
        worksheet = workbook.add_worksheet()

        worksheet.write(0, 0, "id")  # label first row
        worksheet.write(0, 1, "Accession")
        worksheet.write(0, 2, "ReleaseDate")
        row = 1
        csv_reader = csv.reader(f)
        for sample in csv_reader:
            print(sample[0])
            if str(sample[0].strip()) in samples:
                worksheet.write(row, 0, row)  # id
                # worksheet.write(row, 1, sample.strip()) #accession num
                worksheet.write(row, 1, sample[0])  # accession num
                worksheet.write(row, 2, sample[1])  # release date
                # acc = dict2[sample.strip()].pop()
                # worksheet.write(row, 2, getFromAcc(str(acc))) #release date
                row += 1
        workbook.close()
    df = pd.read_excel(wbk)
    csv_file_path = acc_path + "/" + str(part) + "accessionList.csv"
    df.to_csv(csv_file_path, index=False)


# transfers relevant files from a big folder into the correct folders pt0, pt1, pt2
def copyFiles():
    # change range number depending on number of partitions. Eg. 5 partitions, then change to range(5)
    for count in range(num_partitions):
        fasta_dst = args.source + "/seqtrack/pt" + str(count) + "/fasta"
        vcf_dst = args.source + "/seqtrack/pt" + str(count) + "/vcf"
        os.makedirs(fasta_dst)
        os.makedirs(vcf_dst)
        acc_path = (
            args.source
            + "/seqtrack/pt"
            + str(count)
            + "/"
            + str(count)
            + "accessionList.csv"
        )

        with open(acc_path, "r") as f:
            reader = csv.reader(f)
            for row in reader:
                # print(row)  # nid to skip first iteration
                if row[1][0] == "A":
                    continue
                # break
                sample = row[1]
                sample = sample.strip()

                # trf fasta files from src to dst
                src = fastavcfloc + "/fasta/" + str(sample) + ".fasta"
                dst = fasta_dst
                shutil.copy(src, dst)
                # trf vcf files from src1 to dst1
                src1 = fastavcfloc + "/vcf/" + str(sample) + ".vcf"
                dst1 = vcf_dst
                shutil.copy(src1, dst1)


# converts ERR1234567xxx to ERR1234567.lofreq.fasta or ERR1234567.lofreq.vcf, as required by R script
def renameFiles():
    for count in range(num_partitions):
        fastapath = args.source + "/seqtrack/pt" + str(count) + "/fasta/"
        for file in os.listdir(fastapath):
            sample = file[:-6]
            os.rename(fastapath + file, fastapath + sample + ".lofreq.fasta")
        vcfpath = args.source + "/seqtrack/pt" + str(count) + "/vcf/"
        for file in os.listdir(vcfpath):
            sample = file[:-4]
            os.rename(vcfpath + file, vcfpath + sample + ".lofreq.vcf")


os.makedirs(args.source + "/seqtrack", exist_ok=True)
for i in range(num_partitions):
    path = args.source + "/seqtrack/pt" + str(i)
    os.makedirs(path, exist_ok=True)
    makeExcel(path, i)
copyFiles()
renameFiles()

# generate output file for smk
# open("res/5filepartition.done", "x")
