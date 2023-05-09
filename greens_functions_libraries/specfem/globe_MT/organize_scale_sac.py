import sys
import glob
import os

from obspy import read

''' 
    General script to take .sac files from SPECFEM, rename and reorganize them, and scale them for use in MTUQ.   
'''

#Needs to 1) Make new name based on with run directory file is in, 2) Move files to new directory, 3) Scale files
def file_rename(old_name,mij):
    # Goes from {net}.{sta}.{comp}.sem.sac
    # To {net}.{sta}..[RTZ].{Mij}.sac

    dum = old_name.split('/')[-1].split('.')
    net = dum[0]
    sta = dum[1]
    comp = dum[2][-1]
    mij = str(mij)

    new_name = '%s.%s..%s.%s.sac' % (net,sta,comp,mij)
    return(new_name)

def move(eid,run_num,new_path):
    # Moves files from standard simultaneous runs SPECFEM directory structure to MTUQ-compatible directory structure
    # Hardcoded paths for now
    
    run_dir = 'run000' + str(run_num)
    old_path = './specfemglobe_workdir/%s/OUTPUT_FILES/' % run_dir
    
    temp = os.listdir('./specfemglobe_workdir/%s/'%run_dir)
    Mij = temp[any('*.txt' in file for file in temp)].split('.txt')[0]
    Mij = Mij.capitalize()
    
    sac_list = glob.glob(old_path+'*.sac')
    for file in sac_list:
        old_file = file.split('/')[4]
        new_file = file_rename(old_file,Mij)
        print('Moving %s to %s' % (file,new_path+new_file))
        os.rename(file,new_path+new_file)
    
def scale_gf(sc,event):
    # Scales GFs
    sc = float(sc)
    st = read(event)
    print('Scaling %s by %s' % (event,sc))
    st[0].data = st[0].data / sc
    st.write(ev,format='SAC')
    

# Main work starts here
runs = [1,2,3,4,5,6]
eid = 20090407201255351 #could be read from CMTSOLUTION

new_path = './GFs/%s/' % (str(eid))
mk_new_dir = 'mkdir ' + new_path
os.system(mk_new_dir)

for run in runs:
    move(eid,run,new_path)
    

ev_list = glob.glob(new_path + '*.sac')
for ev in ev_list:
    scale_gf(1e15,ev)