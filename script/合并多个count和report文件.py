import os 

# counts files
file_name = [file for file in os.listdir('./') if file.endswith('count')]
counts = 0
for file in file_name:
    # print(file)
    with open(file, 'r') as fr, open('../counts.MM', 'a') as fw:
        fw.write(fr.read())
    counts += 1    
    percent = (counts / len(file_name) ) * 100
    print('\r' + "#" * round(percent) + ">>> %d%% (%d/%d)" %(percent,counts,len(file_name)),end="")
    
# report files
report_files =  [file for file in os.listdir('./') if file.endswith('report')]
counts = 0
for file in file_name:
    # print(file)
    with open(file, 'r') as fr, open('../report.MM', 'a') as fw:
        fw.write(fr.read())
    counts += 1    
    percent = (counts / len(file_name) ) * 100
    print('\r' + "#" * round(percent) + ">>> %d%% (%d/%d)" %(percent,counts,len(file_name)),end="")
print('finish ~') 

  