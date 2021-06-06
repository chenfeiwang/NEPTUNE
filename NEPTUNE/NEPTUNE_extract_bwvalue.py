#!/usr/bin/env python
# author: wangwen
import os, sys
import subprocess

USAGE = 'USAGE: %s <outputfile> <bed> <datapoints> <bw>+' %sys.argv[0]

def generate_subscript():
    """
    subprocess to get signal from bigwig file, direction not considered now
    """

    if not os.path.isfile('tem/tem.py'):

        SubScriptsFHD = open('tem/tem.py', 'w')
        SubScriptsFHD.write("import os, sys\n")
        SubScriptsFHD.write("\n")
        SubScriptsFHD.write("USAGE = '%ss <name> <bedfile> <datapoints> <bigwigfiles>' %ssys.argv[0]\n" %("%", "%"))
        SubScriptsFHD.write("\n")
        SubScriptsFHD.write("datapoints = int(sys.argv[3])\n")
        SubScriptsFHD.write("bwfhds = [x for x in sys.argv[4:]]\n")
        SubScriptsFHD.write("\n")
        SubScriptsFHD.write("rfhds = [open(sys.argv[1]+'_siteprof'+str(i+1), 'w') for i in range(len(bwfhds))]\n")
        SubScriptsFHD.write("bed = {}\n")
        SubScriptsFHD.write("bedfhd = open(sys.argv[2])\n")
        SubScriptsFHD.write("for line in bedfhd.xreadlines():\n")
        SubScriptsFHD.write("    #for i in range(len(bwfhds)):\n")
        SubScriptsFHD.write("    #    rfhds[i].write(line.strip()+'\\t')\n")
        SubScriptsFHD.write("    line = line.strip().split()\n")
        SubScriptsFHD.write("    chrom = line[0]\n")
        SubScriptsFHD.write("    start, end = int(line[1]), int(line[2])\n")
        # SubScriptsFHD.write("    featurename = line[2]\n")
        SubScriptsFHD.write("    for i in range(len(bwfhds)):\n")
        SubScriptsFHD.write("        output = os.popen('bigWigSummary %ss %ss %sd %sd %sd' %s(bwfhds[i], chrom, max(0, start), end, datapoints)).read().strip().split()\n" %("%", "%", "%", "%", "%", "%"))
        SubScriptsFHD.write("        if len(output) == datapoints:\n")
        SubScriptsFHD.write("            output.insert(0,str(end))\n")
        SubScriptsFHD.write("            output.insert(0,str(start))\n")
        SubScriptsFHD.write("            output.insert(0,chrom)\n")
        SubScriptsFHD.write("            rfhds[i].write('\\t'.join(['NA' if value == 'n/a' else value for value in output]) + '\\n')\n")
        SubScriptsFHD.write("        else:\n")
        SubScriptsFHD.write("            output=[]\n")
        SubScriptsFHD.write("            for x in range(datapoints):\n")
        SubScriptsFHD.write("               output.insert(0,'NA')\n")
        SubScriptsFHD.write("            output.insert(0,str(end))\n")
        SubScriptsFHD.write("            output.insert(0,str(start))\n")
        SubScriptsFHD.write("            output.insert(0,str(chrom))\n")
        SubScriptsFHD.write("            rfhds[i].write('\\t'.join([value for value in output]) + '\\n')\n")
        SubScriptsFHD.write("\n")
        SubScriptsFHD.write("bedfhd.close()\n")
        SubScriptsFHD.write("for rfhd in rfhds:\n")
        SubScriptsFHD.write("    rfhd.close()\n")
        SubScriptsFHD.write("\n")
        SubScriptsFHD.close()

def main():

    process = 20
    datapoints = int(sys.argv[3])

    # load bed file
    bed = []
    fhd = open(sys.argv[2])
    for line in fhd.xreadlines():
        bed.append(line)
    fhd.close()

    # creat tem dirctory
    os.system('mkdir tem')

    # add tem files to dirctory, scripts, bwfiles & bed files
    generate_subscript()

    CMD = "cd tem && "
    for i in range(len(sys.argv[4:])):
        CMD += "ln -s ../%s .\n" %(sys.argv[4+i])
    os.system(CMD)

    # tem bed file
    sub_sizes = len(bed) / process
    sub_index = []
    for i in range(process-1):
        sub_index.append((i*sub_sizes, i*sub_sizes+sub_sizes))
    sub_index.append((process*sub_sizes-sub_sizes,len(bed)))
    for i in range(process):
        tembedfhd = open('tem/tem%d_%s.bed' %((i+1), sys.argv[1]), 'w')
        tembedfhd.write("".join(bed[sub_index[i][0]:sub_index[i][1]]))
        tembedfhd.close()

    processes = []
    for i in range(process):
        CMD = "cd tem && python2.7 tem.py %s_tem_%d tem%d_%s.bed %d %s" %(sys.argv[1], i+1, i+1, sys.argv[1], datapoints, ' '.join(sys.argv[4:]))
        # print("Sub Command: ", CMD)
        processes.append(subprocess.Popen(CMD, shell=True))
    stats = []
    for i in range(process):
        stats.append(os.waitpid(processes[i].pid,0))
    CMD = ''
    for i in range(len(sys.argv[4:])):
        CMD += "cat "
        for j in range(process):
            CMD += "tem/%s_tem_%d_siteprof%d " %(sys.argv[1], j+1, i+1)
        CMD += "| gzip - > %s_siteprof%d.gz && \\\n" %(sys.argv[1], i+1)
    CMD = CMD[:-5]
    rfhd = open('tem.sh', 'w')
    rfhd.write(CMD)
    rfhd.close()
    # print("Cat Command: ", CMD)
    p = subprocess.Popen('bash tem.sh', shell = True)
    try:
        os.waitpid(p.pid, 0)
    except OSError:
        pass
    p = subprocess.Popen('rm -r tem tem.sh', shell = True)
    try:
        os.waitpid(p.pid, 0)
    except OSError:
        pass


# program running
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        Info("User interrupts me! ;-) See you!")
        sys.exit(0)
