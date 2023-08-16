#Code for Running Specfem3D for the 6 elementary sources and then, processing the output files for being read by MTUQ. 
#Felix Rodriguez Cardozo, September 2022
#Please refer to the manual before using this algorithm. 


import argparse
from obspy.core import read
from obspy import read
from obspy.geodetics.base import gps2dist_azimuth
from obspy.signal.rotate import rotate_ne_rt
import glob
import os
import numpy as np

def check_solver(id,depth,dir):
    if os.path.exists('{}/OUTPUT_FILES/DATABASES_MPI/'.format(dir)):
        if os.path.exists('SOLVER_REPO/{}/{}'.format(id,depth)):
            list_MT = glob.glob('SOLVER_REPO/{}/{}/M*'.format(id,depth))
            if len(list_MT) == 6:
                return(1)
            else:
                print('Elementary sources directories in SOLVER_REPO/{}/{}/ are not complete. You need 6: Mpp,Mrp,Mrr,Mrt,Mtp,Mtt'.format(id,depth))
                return(0)
        else:
            print('I cannot find the path SOLVER_REPO/{}/{}'.format(id,depth))
            return(0)
    else:
        print('I cannot find the {}/OUTPUT_FILES/DATABASES_MPI/ database'.format(dir))
        return(0)

def check_synthetics(id,depth,MT):
    if os.path.exists('SOLVER_REPO/{}/eventinfo.txt'.format(id)):
        if os.path.exists('SOLVER_REPO/{}/STATIONS'.format(id)):
            cont = 0
            for m in MT:
                check_syn = glob.glob('SOLVER_REPO/{}/{}/{}/OUTPUT_FILES/*sem*'.format(id,depth,m))
                if len(check_syn) == 0:
                    print(' I cannot find synthetics in SOLVER_REPO/{}/{}/{}/OUTPUT_FILES/'.format(id,depth,m))
                    print('It seems you have not run the solver for {} elementary source'.format(m))
                else:
                    cont = cont+1
            if cont == 6:
                return(1)
            else:
                print('Your elementary source simulations are not complete')
                return(0)
        else:
            print('I cannot find the file SOLVER_REPO/{}/STATIONS'.format(id))
            print('Put there the STATIONS file used in Specmfem3D in Geographical Coordinates')
            return(0)
    else:
        print('I cannot find the  file SOLVER_REPO/{}/eventinfo.txt'.format(id))
        print('Make a file like this: YEAR-MONTH-DAYTHOUR:MINUTE:SECONDS LATITUDE LONGITUDE DEPTH MAGNITUDE.\ For example:')
        print('2017-12-01T02:32:44.0 30.7340 57.3900 15.0 5.9')
        return(0)


def run_solvers(id,depth,MT,dir):

    os.chdir('SOLVER_REPO/{}/{}/'.format(id,depth))
    
    for m in MT:
        print('###Running Specfem3D solver for SOLVER_REPO/{}/{}/{}###'.format(id,depth,m))
        os.chdir('{}/OUTPUT_FILES'.format(m))
        s_link = 'ln -s ../../../../../{}/OUTPUT_FILES/DATABASES_MPI/ .'.format(dir)
        print(s_link)
        os.system(s_link)
        os.chdir('../')
        schedule_run = 'sbatch schedule_solver.sh'
        print(schedule_run)
        #os.system(schedule_run)
        os.chdir('../')
        
def process_elementary_sources(MT,id,depth):
    os.chdir('SOLVER_REPO/{}'.format(id))
    print('******Getting event and stations information*******')
    evla,evlo,evdp,time = get_event_info()
    st_info = get_stations_info()
    os.chdir('{}'.format(depth))
    print('###### Converting Specfem3D files into SAC and mini-seed files #####')
    mkdir_output = 'mkdir PROCESSED'
    os.system(mkdir_output)

    for m in MT:
        print('**********Converting to SAC files {} synthetics********'.format(m))
        type_sim,sc = txt2sac(m,time)
        print('********Completing Header Information for {} SAC files ****'.format(m))
        list_sac = glob.glob('PROCESSED/*.sac')
        for s in list_sac:
            stla,stlo = index_st(st_info,s)
            complete_header(s,evla,evlo,evdp,stla,stlo)
        print('********Rotating to R and T {} SAC files ****'.format(m))    
        rotate(m)

    if type_sim == 'semv':
        print('*******Integrating Green Functions*****')
        integrate_gf()

    print('*******Scaling Green Functions*****') 
    scale_gf(sc)

def get_event_info():
    open_info = open('eventinfo.txt','r')
    read_info = open_info.readlines()
    time = read_info[0].split()[0]
    evla = read_info[0].split()[1]
    evlo = read_info[0].split()[2]
    evdp = read_info[0].split()[3]
    return(evla,evlo,evdp,time)

def get_stations_info():
	open_STATIONS = open('STATIONS','r')
	read_STATIONS = open_STATIONS.readlines()
	open_STATIONS.close()
	return(read_STATIONS)

def txt2sac(m,time):

    DIR = '{}/OUTPUT_FILES/*sem*'.format(m)

    list_txt=glob.glob(DIR)

    #Extracting the tail of every sem file for detecting the type of simulation
    type_sim = list_txt[0].split('/')[-1].split('.')[-1]

    if type_sim == 'semv':
        print('********** I detected velocity simulations **********')
    elif type_sim == 'semd':
        print('********** I detected displacement simulations **********')

    #Read CMTSOLUTION file for guessing the scaling factor
    open_cmtsolution= open('{}/DATA/CMTSOLUTION'.format(m),'r')
    read_cmtsolution = open_cmtsolution.readlines()
    
    for line in read_cmtsolution:
        aux = line.split()[0].split(':')[0]
        if aux == m:
            sa = float(line.split()[-1])
            sc = sa/1e+7

    for i in range(len(list_txt)):
        print('Pre-processing {}'.format(list_txt[i]))
        open_data=open(list_txt[i],'r')
        data=open_data.readlines()
        open_data.close()

        samples = len(data)
        time_s = np.round(np.abs(float(data[0].split()[0])) + np.abs(float(data[-1].split()[0])))
        sps = np.round(float(samples/time_s))
   
        data_series=[]
        exp=[]
        for j in range(len(data)):
            data_series.append(data[j].split()[1])
            #exp.append(data_series[-1].split('E-')[1])
            #print data_series
        write_data=open(list_txt[i]+'.ascii','w')
        name_stat=list_txt[i].split('/')[-1]
        split_name=name_stat.split('.')
        name_st=split_name[1]
        component=split_name[2]
        network=split_name[0]
        name_stat=name_st+'_'+component+'_00'+'_'+network+'_R'
        #header='TIMESERIES '+name_stat+','+' 42000 samples, 100 sps, 2017-12-01T02:32:44.000000, SLIST, FLOAT,\n'
        header='TIMESERIES '+name_stat+','+' {} samples, {} sps, {}, SLIST, FLOAT,\n'.format(samples,sps,time)
        write_data.write(header)
        for j in range(len(data_series)):
            #aux_data=float(data_series[j])*float('1e'+exp[j])
            write_data.write(data_series[j]+'\n')
        write_data.close()
        st = read(list_txt[i]+'.ascii')

        new_name = rename(list_txt[i],m)
        st.write(new_name, format='SAC')
        print (st)

        remove_ascii ='rm {}.ascii'.format(list_txt[i]) 
        print(remove_ascii)
        os.system(remove_ascii)

        mv_sac = 'mv {} PROCESSED/'.format(new_name)
        print(mv_sac)
        os.system(mv_sac)
        print('\n')
    
    return(type_sim,sc)

def rename(old_name,m):
    #Mtt/OUTPUT_FILES/IR.JHBN.HXZ.semv
    #Mtt/OUTPUT_FILES/IR.CHBR..Z.Mrp.sac
    aux=old_name.split('/')[-1].split('.')
    network = aux[0]
    st = aux[1]
    comp = aux[2][-1]
    if comp == 'X':
        comp = 'E'
    elif comp == 'Y':
        comp = 'N'
    new_name = '{}/OUTPUT_FILES/{}.{}..{}.{}.sac'.format(m,network,st,comp,m)
    return(new_name)

def index_st(st_info,s):
    aux = s.split('/')[-1]
    st_name = aux.split('.')[1]
    for i in st_info:
        if st_name in i :
            pos = st_info.index(i)
            break
    stla = st_info[pos].split()[2]
    stlo = st_info[pos].split()[3]
    return(stla,stlo)

def complete_header(s,evla,evlo,evdp,stla,stlo):
    #PROCESSED/IR.CHBR..Z.Mrp.sac
    st = read(s)
    aux = s.split('/')[-1]
    network = aux.split('.')[0]

    st[0].stats.sac.evla = evla
    st[0].stats.sac.evlo = evlo
    st[0].stats.sac.evdp = evdp
    st[0].stats.sac.stla = stla
    st[0].stats.sac.stlo = stlo
    st[0].stats.sac.knetwk = network
    st[0].stats.sac.khole = ''
    st[0].stats.network = network
    st[0].stats.location = ''
    st[0].stats.lcalda ='TRUE'

    st.write(s,format='SAC')

def rotate(m):
	#PROCESSED/IR.CHBR..Z.Mrp.sac
	list_st = glob.glob('PROCESSED/*..Z.*.sac'.format(m))
	for t in list_st:
		stZ = read(t)
		net = stZ[0].stats.network
		source_lat = stZ[0].stats.sac.stla
		source_lon = stZ[0].stats.sac.stlo
		st_lat = stZ[0].stats.sac.evla
		st_lon = stZ[0].stats.sac.evlo
		baz = gps2dist_azimuth(source_lat,source_lon,st_lat,st_lon)
        
		st_name = t.split('.')[1]
		print('Rotating for {},{}'.format(m,st_name))
		stN = read('PROCESSED/{}.{}..N.{}.sac'.format(net,st_name,m))
		stE = read('PROCESSED/{}.{}..E.{}.sac'.format(net,st_name,m))

		rotated = rotate_ne_rt(stN[0].data,stE[0].data,baz[1])
    
		stN[0].data = rotated[0]
		stN[0].stats.sac.kcmpnm = 'BHR'
		stN[0].stats.channel = 'BHR'
		stN[0].stats.sac.cmpaz = baz[2]

		stE[0].data = rotated[1]
		stE[0].stats.sac.kcmpnm = 'BHT'
		stE[0].stats.channel = 'BHT'

		if baz[2] + 90 >= 360:
			stE[0].stats.sac.cmpaz = baz[2] + 90 - 360
		else :
			stE[0].stats.sac.cmpaz = baz[2] + 90

		r_name = 'PROCESSED/{}.{}..R.{}.sac'.format(net,st_name,m)
		t_name = 'PROCESSED/{}.{}..T.{}.sac'.format(net,st_name,m)
        
		stN.write(r_name,format='SAC')
		stE.write(t_name,format='SAC')

def integrate_gf():
    list_ev = glob.glob('PROCESSED/*.sac')

    for ev in list_ev:
        print('Integrating {}'.format(ev))
        st = read(ev)
        st[0].data = np.cumsum(st[0].data)*st[0].stats.delta
        st.write(ev,format='SAC')

def scale_gf(sc):
    list_ev = glob.glob('PROCESSED/*.sac')

    for ev in list_ev:
        print('Scaling {}'.format(ev))
        st = read(ev)
        st[0].data = st[0].data/sc
        st.write(ev,format='SAC')

parser = argparse.ArgumentParser()
parser.add_argument("-o",type=str,help='Options: s (run solver), p (process output files). E.g., python make_specfem3D_GF.py  -o s, will run the Specfem3D solver for the 6 elementary sources. ')
parser.add_argument("-ev",type=str,help="Event ID. E.g. -ev  20171201023244")
parser.add_argument("-ed",type=str,help="Earthquake depth in km E.g. -ed  15")
parser.add_argument("-dir",type=str,help="Database dir E.g. -dir DATABASES_REPO/1D_VM/16_PROC/")

args = parser.parse_args()
opt= str(args.o)
id = str(args.ev)
depth = str(args.ed)
dir = str(args.dir)

cwd = os.getcwd()
MT = ['Mpp','Mrp','Mrr','Mrt','Mtp','Mtt']
go = check_solver(id,depth,dir)

if go == 0:
    print('End')
else:
    if opt == 's':
        run_solvers(id,depth,MT,dir)
        os.chdir(cwd)
    elif opt == 'p':
        print('Pre-processing Specfem3D output files')
        go2 = check_synthetics(id,depth,MT)
        if go2 == 0:
            print('End')
        else:
            process_elementary_sources(MT,id,depth)
    else: 
        print('Option {} was not recognized'.format(opt))