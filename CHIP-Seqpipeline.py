#!/usr/bin/python

_version_ = '1.0'

import os
import sys
import argparse
import commands
import time
import shutil

def get_args():
	tool = os.path.basename(sys.argv[0])
	parser = argparse.ArgumentParser(description='programm: '+ tool + ' (CHIP-Seq data analyse pipeline) \n\nversion: '+ _version_+'\nAuthor:\tzp' + '\nEmail:zhangpeng1334880@gmail.com',prog=tool,formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument('-v','--version',action='version',version='%(prog)s v'+_version_)

	parser.add_argument('-alignment_type',dest='alignment_type',help='align the reads to reference genome.(default:bwa mem)\nrecomand(bwa-mem,bowtie,bowtie2)',default='bwa mem')
	parser.add_argument('-alignment_index',dest='alignment_index',help='corresponding index file\n(keep constistent with -alignment_type eg:if your choose bwa this is bwa index)',default='none')
	parser.add_argument('-O',dest='output',help='The output directory to put results and log files in',default='OUTPUT')
	parser.add_argument('reads1',help='Files with #1 mates, paired with files in <m2>.\nCould be gzip\'ed (extension: .gz) or bzip2\'ed (extension: .bz2).',default='none')
	parser.add_argument('reads2',help='Files with #2 mates, paired with files in <m1>.\nCould be gzip\'ed (extension: .gz) or bzip2\'ed (extension: .bz2).',default='none')

	args = parser.parse_args()
	if args.alignment_type == 'bwa-mem':
		args.alignment_type = 'bwa mem' 
	return args

def check_args():
	if args.alignment_index == 'none':
		print '\tWARNING! No parameter {0}'.format('-alignment_index') 
		return False
	if len(args.reads1) != len(args.reads2):
		print '\tWARNING! reads count no not match!'
		exit()
	if len(args.reads1) ==0 or len(args.reads2) == 0:
		print '\tWARNING! no reads input, please check {0} and {1}'.format('reads1','reads2')
		exit()
	return True

def create_dir(dir):
    if os.path.isdir(dir):
    	if os.listdir(dir):
    		print 'WARNING directory %s is not empty, please change a directory!' % dir
    	exit()
    	# shutil.rmtree(dir)
    os.mkdir(dir)

def alignment(sample_name,file_in,file_out):

	command = args.alignment_type + ' ' + args.alignment_index + ' ' + file_in + ' ' + file_out + ' > ' + args.output +'/' + sample_name + '.sam'
	text = 'alignment {0} {1} and {2} is over'.format(args.alignment_type,reads1[i],reads2[i]) 
	print text + ' '*(width-len(text)) + time.strftime('%X',time.localtime())
	# print 'bwa mem {0} and {1} is over		 {2} '.format(file_in,file_out,time.strftime('%X',time.localtime()) )
	commands.getstatusoutput(command)
	return True

def sam_to_bam(file_in,file_out):
	command = 'samtools view -bS {0} > {1}/{2}'.format(file_in,args.output,file_out)
	commands.getstatusoutput(command)
	return True

def sam_sort(file_in,file_out):
	command = 'samtools sort {0} {1}{2}'.format(file_in,args.output,file_out)
	commands.getstatusoutput(command)
	return True

def samtools_index(file_in):
	command = 'samtools index {0}'.format(file_in)
	commands.getstatusoutput(command)
	return True

def samtools_flagstat(file_in):
	command = 'samtools flagstat {0}'.format(file_in)
	output = commands.getoutput(command)
	# f = open('flagstat.txt','w')
	# f.write(output)
	return True

def mark_duplicates(file_in):
	command = 'java -jar MarkDuplicates.jar INPUT={0}.sorted.bam OUTPUT={1}/{0}.rmdup.bam METRICS_FILE={0}.metrics'.format(file_in,args.output)
	commands.getstatusoutput(command)
	return True

def premacs2(sample_name):
	name = sample_name
	if sam_to_bam(name + '.sam',name + '.bam') == True:
		text = 'sam to bam {0} and {1} is over'.format(reads1[i],reads2[i]) 
		print text + ' '*(width-len(text)) + time.strftime('%X',time.localtime())

	if sam_sort(name + '.bam',name + '.sorted' ) == True:
		text = 'samtools sort {0} and {1} is over'.format(reads1[i],reads2[i]) 
		print text + ' '*(width-len(text)) + time.strftime('%X',time.localtime())

	if samtools_index(name + '.sorted.bam' ) == True:
		text = 'samtools index {0} and {1} is over'.format(reads1[i],reads2[i]) 
		print text + ' '*(width-len(text)) + time.strftime('%X',time.localtime())

	if samtools_flagstat(name + '.sorted.bam') == True:
		text = 'samtools flagstat {0} and {1} is over'.format(reads1[i],reads2[i]) 
		print text + ' '*(width-len(text)) + time.strftime('%X',time.localtime())

	if mark_duplicates(name):
		text = 'picard mark_duplicates {0} and {1} is over'.format(reads1[i],reads2[i]) 
		print text + ' '*(width-len(text)) + time.strftime('%X',time.localtime())

	return True

def run_macs2(treat,control='none',name='none'):
	treat = treat.split(',')
	if len(treat) == 1:
		command = 'macs2 callpeak -g hs -f BAM -t {0}.rmdup.bam -c {1}.rmdup.bam -n {2} --output_dir {2}_output/'.format(treat,control,name)	
		commands.getstatusoutput(command)	
		text = 'macs2 callpeak {0} and {1} is over'.format(str(treat[0]),control) 

	if len(treat) == 4 and control == 'none':
		command = 'macs2 callpeak -g hs -f BAM -t {0}.rmdup.bam {1}.rmdup.bam {2}.rmdup.bam {3}.rmdup.bam -n {4} --output_dir {4}_output/'.format(treat[0],treat[1],treat[2],treat[3],name)
		commands.getstatusoutput(command)
		text = 'macs2 callpeak {0} is over'.format(','.join(treat)) 

	print text + ' '*(width-len(text)) + time.strftime('%X',time.localtime())
	return True

def peak(file_a,file_b,peak_out):
	name = 'peak_condition_' + file_a + file_b
	command = 'bedtools intersect -wa {0} -wb {1} > {2}/{3}'.format(file_a,file_b,args.output,peak_out)
	commands.getstatusoutput(command)
	text = 'peak {0} and {1} is over'.format(file_a,file_b) 
	print text + ' '*(width-len(text)) + time.strftime('%X',time.localtime())
	return True

if __name__ == '__main__':

	args = get_args()
	flag = check_args()
	if flag == True:

		reads1 = args.reads1.split(',')
		reads2 = args.reads2.split(',')
		width = 62
		create_dir(args.output)

		print '*'*100
		print 'STEP1:	alignment\n\nalignment type: ' + args.alignment_type

		for i in range(len(reads1)):
			sample_name = reads1[i].split('_')[0]
			if alignment(sample_name,reads1[i],reads2[i]) != True:
				print 'STOP! alignment failed, please check the alignment process!'


		print '*'*100
		print 'STEP2:	pre-MACS\n'
		sample_names = []
		for i in range(len(reads1)):
			sample_name = reads1[i].split('_')[0]
			sample_names.append(sample_name)

			if premacs2(sample_name) == True:
				pass
			else:
				print 'STOP! pre-MACS failed,please check the pre-MACS process!'


		print '*'*100
		print 'STEP3:	MACS call peaks\n'

		treat_control_list = { sample_names[0]:sample_names[1],sample_names[2]:sample_names[3],sample_names[4]:sample_names[5],sample_names[6]:sample_names[7] }
		for treat in treat_control_list.keys():
			if run_macs2(treat,control=treat_control_list[treat],name='peak_' + treat) == True:
				pass
			else:
				print 'STOP! macs2 failed, please check the file {0} and {1}!'.format(treat,treat_control_list[treat])


		print '*'*100
		print 'STEP4:	Peak\n'

		condition1 = [ sample_names[0],sample_names[1],sample_names[2],sample_names[3] ]
		condition2 = [ sample_names[4],sample_names[5],sample_names[6],sample_names[7] ]
		condition = [','.join(condition1),','.join(condition2)]
		for i in range(len(condition)):
			if run_macs2(condition[i],name='peak_all_' + str(i+1)) == True:
				pass
			else:
				print 'STOP! macs2 failed, please check the file {0}!'.format(condition[i])


		print '*'*100
		print 'STEP5:	compare\n'

		name = treat_control_list.keys()
		##condition1
		peak('peak_'+name[0],'peak_all_1','peak_ln_1')
		peak('peak_'+name[1],'peak_all_1','peak_ln_2')
		peak('peak_ln_1','peak_ln_2','peak_condition_1')

		###condition 2
		peak('peak_'+name[2],'peak_all_2','peak_ln_3')
		peak('peak_'+name[3],'peak_all_2','peak_ln_4')
		peak('peak_ln_4','peak_ln_4','peak_condition_2')

		print '*'*100
		text = 'Game over!'
		print text + ' '*(width-len(text)) + time.strftime('%X',time.localtime()) + '\n'
		print '*'*100

	else:
		exit()
# !python CHIP-Seqpipeline.py -alignment_index BWAIndex/ -O OUTPUT treat1_1,input1_1,treat2_1,input2_1,treat3_1,input3_1,treat4_1,input4_1 treat1_2,input1_2,treat2_2,input2_2,treat3_2,input3_2,treat4_2,input4_2

# !python CHIP-Seqpipeline.py -h

