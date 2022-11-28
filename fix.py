import csv

vals = []
with open("endf_ar.txt", "r") as file:
    reader = csv.reader(file, delimiter=",")
    for row in reader:
        temp_energy = ""
        temp_cross_section = ""
        energy = True
        for c in row[0]:
            if(c == " " or c == "\n"):
                if(energy):
                    energy = False
                continue
            else:
                if(energy):
                    temp_energy += c
                else:
                    temp_cross_section += c
        vals.append([float(temp_energy),float(temp_cross_section)])

with open("endf_ar.csv", "w") as file:
    writer = csv.writer(file, delimiter=",")
    writer.writerows(vals)