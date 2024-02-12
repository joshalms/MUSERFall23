fp = 'links_3.txt'
with open(fp, 'r') as f:
    for line in f:
        cf= ','.join(part.strip() for part in line.split(','))