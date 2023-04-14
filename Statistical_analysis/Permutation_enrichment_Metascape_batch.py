import os
from multiprocessing.pool import Pool
input_dir='data/Genelist/WMH_pos/'
output_dir=input_dir.replace('Genelist','gene_result_new')
start_num=25
def mkdir(p):
    if not os.path.exists(p):
        os.makedirs(p)
        os.chmod(p,511)
mkdir(output_dir)
def start(i):
    os.system('bin/up.sh '+str(i))
def stop(i):
    os.system('bin/down.sh '+str(i))
def process(input_list,i):
    #start(i)
    for t in input_list:
        t_out=os.path.join(output_dir,t[:-5])
        if os.path.exists(t_out):
            if len(os.listdir(t_out))>10:
                continue
        cmd_line='bin/ms.sh {t_num} -o /{t_result} /{t} --option /data/example/option.json'.format(t_num=i,t=os.path.join(input_dir,t),t_result=t_out)
        os.system(cmd_line)
    stop(i)
if __name__=='__main__':
    input_list=[]
    d_num=3
    for i in range(d_num):
        start(i+start_num)
    for t in os.listdir(input_dir):
        if t.endswith('.xlsx'):
            input_list.append(t)
    s_num=int(len(input_list)/d_num)
    new_input=[]
    for i in range(d_num):
        if i==d_num-1:
            args=input_list[i*s_num:],i+start_num
        else:
            args=input_list[i*s_num:i*s_num+s_num],i+start_num
        new_input.append(args)
    p=Pool(d_num)
    p.starmap(process,new_input)
    p.close()
    p.join()
