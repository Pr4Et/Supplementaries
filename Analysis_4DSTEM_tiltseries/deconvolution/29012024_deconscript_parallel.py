#!/usr/bin/python3

import os
from time import perf_counter
import fnmatch
import shutil
import subprocess
from multiprocessing import Pool                                                
import multiprocessing
from numpy import *

####################################################################################################
# This script needs a 3d volume, a tilt file (not the .rawtlt but the refined tilt angles(.tlt))
# and a precomputed PSF.mrc (with the right Pixelsize)!
# e.g. the 3d volume is called "thisisthetomogram_rec.mrc", then the tilt file has to be called
# "thisisthetomogram_rec.tlt". For dual_axis to work, the .tlt files have to be called
# "thisisthetomogram_rec_a.tlt" and "thisisthetomogram_rec_b.tlt".
# This script need IMOD and PRIISM to be called
# It needs a file with some input values! Please adapt as needed 
# For dual_axis it also need the patchcorr log file from IMOD to get the real Axis_rotationangle.
# It needs to be called "thisisthetomogram_rec_patchcorr.log"
##################################################################

# this is tailored to running on the Weizmann wexac cluster
def environment():
    clip_pr=shutil.which("clip")
    ftransform_pr=shutil.which("FTransform3D")
    header_pr=shutil.which("header")
    if not str('clip') in str(clip_pr):
        print('####################################')
        print('### IMOD not loaded ... I will exit the script ###')
        print('load priism with "module load cuda/10.2 IMOD/4.12.10"')
        print('and "export LIBGL_ALWAYS_INDIRECT=0"')
        print('####################################')
        exit()
    if not str('IMOD') in str(header_pr): 
        print('#####################################')
        print('### Script needs the IMOD header ... I will exit the script ###')
        print('## load priism first and IMOD last ##')
        print('#####################################')
        exit() 
    if not str('FTransform3D') in str(ftransform_pr):
        print('####################################')
        print('### PRIISM NOT LOADED ... I will exit the script ###')
        print('load priism with "module load priism/4.7.2" and "source /apps/RH7U2/general/priism/4.7.2/Priism_setup.sh"')
        print('####################################')
        exit()
    else: 
        print('####################################')
        print('############# Off I go #############')
        print('####################################')
        

def write_log_file(logfile):
  if logfile in os.listdir('.'): 
    with open(logfile, 'a') as logf:
      logf.write("\n\n\n\n\n\n")
      logf.write("##### I will continoue another run #####\n")
      logf.write("#####  I have used the " + inp['script'] + " #####\n")
  if not logfile in os.listdir('.'): 
    with open(logfile, 'w') as logf:
      logf.write("##### THIS IS THE LOG OF THE DEVONVOLUTION OF '% s' #####\n" % inp['tomoname'] )
      logf.write("#####  in % s #####\n" % inp['workdir'])
      logf.write("#####  I have used the " + inp['script'] + " #####\n")
    
def input_files(file_to_work_on, path):
    scriptname = '29012024_deconscript_parallel.py'
    orifile = os.path.join(path, file_to_work_on)
    tomoname = file_to_work_on.split('.')[0]
    invtomo = tomoname + '_inv.mrc'
    scaletomo = tomoname + '_inv_scaled.mrc'
    logfile = 'decon_' + tomoname + '.log'
    workdir = os.path.join(path, 'deconv_' + tomoname)
    Fan_file = 'FanOut_' + tomoname + '.mrc'
    OTF_file = Fan_file.split('.')[0] + '.otf'
    newfile = os.path.join(workdir, file_to_work_on)
    filenames = {'parentdir': path, 'script': scriptname, 'orifile': orifile, 'tomoname': tomoname, 'invtomo': invtomo, 'scaletomo': scaletomo, 'logfile': logfile, 'workdir': workdir, 'Fan_file': Fan_file, 'OTF_file': OTF_file, 'newfile': newfile}
    if not os.path.exists(workdir):
        os.mkdir(workdir)
        os.chdir(workdir)
    else:
        os.chdir(workdir)
    shutil.copy2(orifile, newfile)
    shutil.copy2(os.path.join(parentdir + '/' + scriptname), workdir + '/' + scriptname)
    shutil.copy2('../input_values.txt', workdir + '/input_values.txt')
    return filenames

def de_para():
    decon_para = {}
    with open(inp['logfile'], 'a') as logf:
        psf_file = os.path.join('/storwis/Labs/cpelbaum/seifer/deconvolution' +'/PSF3x3pipe.mrc') # '/PSF_' + test_dict['pix_nm'] + 'nm.mrc')
        decon_para['PSF'] = psf_file
        with open('input_values.txt', 'r') as inpf:
            for lines in inpf:
                if 'Number_of_Iterations' in lines:
                    iterations = lines.split()[2]
                    ni = iterations.split(',')
                    decon_para['ni'] = ni
                if 'Smoothingfactor' in lines:
                    smf = str(lines.split()[2])
                    smooth = smf.split(',')
                    decon_para['smooth'] = smooth
                if 'Dual_Axis' in lines:
                    if 'Yes' in lines:
                        logf.write("##### I will do a dual_axis deconvolution #####'\n'")
                        a_tlt_file = inp['tomoname'] + '_a.tlt'
                        shutil.copy2(inp['parentdir'] + '/' + a_tlt_file, inp['workdir'] + '/' + a_tlt_file)
                        b_tlt_file = inp['tomoname'] + '_b.tlt'
                        shutil.copy2(inp['parentdir'] + '/' + b_tlt_file, inp['workdir'] + '/' + b_tlt_file)
                        patchcor = inp['tomoname'] + '_patchcorr.log'
                        shutil.copy2(inp['parentdir'] + '/' + patchcor, inp['workdir'] + '/' + patchcor)
                        decon_para['tlt_file'] = a_tlt_file
                        decon_para['b_tlt_file'] = b_tlt_file
                        decon_para['patchcor'] = patchcor
                        decon_para['dual_axis'] = 'Yes'
                    elif 'No' in lines:
                        logf.write("##### I will do a single_axis deconvolution #####'\n'")
                        tlt_file = inp['tomoname'] + '.tlt'
                        decon_para['tlt_file'] = tlt_file
                        shutil.copy2(inp['parentdir'] + '/' + tlt_file, inp['workdir'] + '/' + tlt_file)
                        decon_para['dual_axis'] = 'No'
                    else:
                        print('Please write "Dual_Axis = Yes" or "Dual_Axis = No" in the input_values.txt')
                        exit()
                if 'Auxiliary_smoothing' in lines:
                    aux_sm = str(lines.split()[2])
                    lamsmooth = aux_sm.split(',')
                    decon_para['lamsmooth'] = lamsmooth
                if 'Do_summed_FFT' in lines:
                    if 'Yes' in lines:
                        decon_para['do_fft'] = 'Yes'
                    elif 'No' in lines:
                        decon_para['do_fft'] = 'No'
                    else:
                        print('Please write "Do_summed_FFT = Yes" or "Do_summed_FFT = No" in the input_values.txt')
                        exit()
                if 'Do_summed_PS' in lines:
                    if 'Yes' in lines:
                        decon_para['do_ps'] = 'Yes'
                    elif 'No' in lines:
                        decon_para['do_ps'] = 'No'
                    else:
                        print('Please write "Do_summed_PS = Yes" or "Do_summed_PS = No" in the input_values.txt')
                        exit()
                threads = multiprocessing.cpu_count()
                decon_para['np'] = str(threads)
                if 'Memory_MB' in lines: 
                    memory = str(lines.split()[2])
                    decon_para['mem'] = memory
                if 'Needs_Invertion' in lines: 
                    if 'Yes' in lines: 
                        decon_para['do_inv'] = 'Yes'
                    elif 'No' in lines: 
                        decon_para['do_inv'] = 'No'
                    else:
                        print('Please write "Needs_Invertion = Yes" or "Needs_Invertion = No" in the input_values.txt')
                        exit()
            logf.write('\n### These are the values I read from the input_values.txt: ###\n')
            logf.write(str(decon_para) + '\n\n')
    return decon_para

def extract_data(inputfile):
    with open(inp['logfile'], 'r+') as logf:
        tomoheader = subprocess.run(['header', inp['newfile']], capture_output=True, text=True)  # read the header into stdout
        logf.write(tomoheader.stdout)  # write the header into the logfile
        extracted_values = {}  # create a dictonary
        extracted_values.clear()
        # extract the minimum density value
        mini = subprocess.run(['grep', 'Minimum'], capture_output=True, text=True, input=tomoheader.stdout)
        mini = str(mini.stdout).split()[-1]
        extracted_values['mini'] = mini
        # extract the maximum density value
        maxi = subprocess.run(['grep', 'Maximum'], capture_output=True, text=True, input=tomoheader.stdout)
        maxi = str(maxi.stdout).split()[-1]
        extracted_values['maxi'] = maxi
        # extract the pixel size value
        pix_a = subprocess.run(['grep', 'Pixel spacing'], capture_output=True, text=True, input=tomoheader.stdout)
        pix_a = round(float(pix_a.stdout.split()[-1]),2)
        pix_nm = round((float(pix_a) / 10),3)
        extracted_values['pix_A'] = str(pix_a)
        extracted_values['pix_nm'] = str(pix_nm)
        # extract the Image size
        image_size = subprocess.run(['grep', 'Number of columns'], capture_output=True, text=True, input=tomoheader.stdout)
        xy = str(image_size.stdout).split()[-2]
        z = str(image_size.stdout).split()[-1]
        extracted_values['xy'] = xy
        extracted_values['z'] = z
        # write a summary into the logfile
        logf.write('\n### These are the values I read from the header of the tomogram: ###\n')
        logf.write('Minmum = ' + mini + '\n')
        logf.write('Maximum = ' + maxi + '\n')
        logf.write('Pixel size (A) = ' + str(pix_a) + '\n')
        logf.write('Pixel size (nm) = ' + str(pix_nm) + '\n')
        logf.write('Image size is = ' + xy + ',' + xy + ',' + z + '\n')
    return extracted_values

def run_process(process):                                                             
    subprocess.run(process, shell=True)  
    
def fan_gen(n_cores):
    with open(inp['logfile'], 'a') as logf:
        deg = genfromtxt(decon_para['tlt_file'])
        params = []
        for i in range(0,len(deg)):
            rot_cmd = 'rotatevol -memory ' + decon_para['mem'] + ' -FillValue 0 -OutputSizeXYZ ' + test_dict['xy'] + ','+ test_dict['xy'] + ',' + test_dict['z'] + ' -a 0,' + str(deg[i]) + ',0 -i ' + decon_para['PSF'] + ' -ou rotvol/PSF_' + str(deg[i]) + 'degY.mrc'
            logf.write(rot_cmd + '\n')
            params.append(rot_cmd)
        processes = ()
        for i in range(0,len(params)):
            processes = processes + ((params[i]),)
        pool = Pool(processes = n_cores)                                                        
        pool.map(run_process, processes)       
        pool.close()     
        #Now I add up all the tilted PSFs to the a axis FanOut
        cadd_cmd = 'clip add -3d rotvol/PSF_*.mrc ' + 'rotvol/' + inp['Fan_file']
        subprocess.run(cadd_cmd, shell=True)
        logf.write(cadd_cmd + '\n')

def fan_gen_b(n_cores):
    with open(inp['logfile'], 'a') as logf:
        with open(decon_para['patchcor'], 'r') as patchcor:
            for lines in patchcor:
                if 'AxisRotationAngle' in lines:
                    axisrotationangle = lines.split()[-1]
                    logf.write('The AxisrotationAngle is: ' + axisrotationangle + '\n\n')
                    axisrot_dif = str(round(90 - float(axisrotationangle), 2))
        deg = genfromtxt(decon_para['b_tlt_file'])
        params = []
        for i in range(0,len(deg)):
            rot_cmd = 'rotatevol -memory ' + decon_para['mem'] + ' -FillValue 0 -OutputSizeXYZ ' + test_dict['xy'] + ','+ test_dict['xy'] + ',' + test_dict['z'] + ' -a 0,' + str(deg[i]) + ',0 -i ' + decon_para['PSF'] + ' -ou rotvol_b/PSF_' + str(deg[i]) + 'degY.mrc'
            logf.write(rot_cmd + '\n')
            params.append(rot_cmd)
        processes = ()
        for i in range(0,len(params)):
            processes = processes + ((params[i]),)
        pool = Pool(processes = n_cores)                                                        
        pool.map(run_process, processes)       
        pool.close()
        #Now I add up all the tilted PSFs to the b axis FanOut
        cadd_cmd = 'clip add -3d rotvol_b/PSF_*.mrc ' + 'rotvol_b/' + inp['Fan_file']
        baxis_flipxy_cmd = 'clip flipxy -3d rotvol_b/' + inp['Fan_file'] + ' rotvol_b/' + inp['Fan_file']
        baxis_roty_cmd = 'rotatevol -a 0,' + axisrot_dif + ',0 -i rotvol_b/' + inp['Fan_file'] + ' -ou rotvol_b/' + inp['Fan_file']
        cmds = [cadd_cmd, baxis_flipxy_cmd, baxis_roty_cmd]
        for cmd in cmds:
            subprocess.run(cmd, shell=True)
            logf.write(cmd + '\n')


def fan_prep():
    with open(inp['logfile'], 'a') as logf:
        if not os.path.isfile(inp['Fan_file']):
            if decon_para['dual_axis'] == 'Yes':
                logf.write("##### I will create a new dual_axis FanOut_PSF with these commands: #####'\n'")
                if not os.path.exists('rotvol'):
                    os.mkdir('rotvol')
                if not os.path.exists('rotvol_b'):
                    os.mkdir('rotvol_b')
                fan_gen(int(decon_para['np']))
                fan_gen_b(int(decon_para['np']))
                cadd_ab_cmd = 'clip add -3d rotvol/' + inp['Fan_file'] + ' rotvol_b/' + inp['Fan_file'] + ' ' + inp['Fan_file']
                subprocess.run(cadd_ab_cmd, shell=True)
                logf.write('I have made a dual_axis PSF!!!' + '\n')
                logf.write(cadd_ab_cmd + '\n\n')
            if decon_para['dual_axis'] == 'No':
                logf.write("##### I will create a new single_axis FanOut_PSF with these commands: #####'\n'")
                if not os.path.exists('rotvol'):
                    os.mkdir('rotvol')
                fan_gen(int(decon_para['np']))
                logf.write('I have made a single_axis PSF!!!' + '\n')
                shutil.copy2('rotvol/' + inp['Fan_file'], inp['Fan_file'])
            #flipz_cmd = 'clip flipz ' + inp['Fan_file'] + ' ' + inp['Fan_file']
            #subprocess.run(flipz_cmd, shell=True)
            #logf.write(flipz_cmd + '\n\n')
            scale_cmd = 'newstack -float 2 -meansd ' + scale + inp['Fan_file'] + ' -ou ' + inp['Fan_file']
            subprocess.run(scale_cmd, shell=True)
            logf.write("##### I have scaled the density and changed the mode to 32bit floating point with this command: #####'\n'")
            logf.write(scale_cmd + '\n\n')
            alt_head_cmd = 'alterheader -RootMeanSquare ' + inp['Fan_file']
            subprocess.run(alt_head_cmd, shell=True, stdout=logf)
        else:
            logf.write('I will use the existing FanOut file: ' + inp['Fan_file'] + '\n\n')
        logf.write("##### The header of the FanOut #####'\n'")
        subprocess.run(['header', inp['Fan_file']], stdout=logf)

def ftransform():
    with open(inp['logfile'], 'a') as logf:
        if not os.path.isfile(inp['OTF_file']):
            # COMMAND TO RUN FTransform 3D ####
            new_xy = str(int(test_dict['xy']) - 1)
            new_z = str(int(test_dict['z']) - 1)
            if int(test_dict['xy']) % 2 == 0:
                shift_xy = int(test_dict['xy']) / 2
            else:
                shift_xy = (int(test_dict['xy']) - 1) / 2
            shift_xy = str(round(shift_xy, 0))
            if int(test_dict['z']) % 2 == 0:
                shift_z = int(test_dict['z']) / 2
            else:
                shift_z = (int(test_dict['z']) - 1) / 2
            shift_z = str(round((shift_z), 0))
            shift = shift_xy + ':' + shift_xy + ':' + shift_z

            logf.write("##### I have created a OFT with this command: #####'\n'")
            # this command writes a file called OTF_file
            ft_cmd = 'FTransform3D ' + inp['Fan_file'] + ' ' + inp['OTF_file'] + ' -x1=0:' + new_xy + ' -y1=0:' + new_xy + ' -z1=0:' + new_z + ':1 -w1=0 -t1=0:0:1 -res1=0 -mode1=complex -title='' -real_complex -same_units -shift=' + shift + ' -xpad=default -xpad_value=0:0 -ypad=default -ypad_value=0:0 -zpad=default -zpad_value=0:0 -3d=z'
            logf.write(ft_cmd + '\n\n')
            subprocess.run(ft_cmd, shell=True)
        else:
            logf.write('I will use the existing OTF file: ' + inp['OTF_file'])
    
def invert_input_file(inputfile):
    with open(inp['logfile'], 'a') as logf:
        inv_cmd = 'newstack -scale ' + test_dict['maxi'] + ',' + test_dict['mini'] + ' -in ' + inputfile + ' -ou ' + inp['invtomo']  # Command to invert the density
        logf.write("##### I have inverted the density with this command: #####'\n'")
        logf.write(inv_cmd + '\n\n')
        subprocess.run(inv_cmd, shell=True)
        logf.write("##### The header of the inverted tomogram #####'\n'")
        subprocess.run(['header', inp['invtomo']], stdout=logf)


def prepare_input_files(inputfile):
    with open(inp['logfile'], 'a') as logf:
        if not os.path.isfile(inp['scaletomo']):
            if decon_para['do_inv'] == 'Yes':
                invert_input_file(inputfile)
            elif decon_para['do_inv'] == 'No':
                inp['invtomo'] = inputfile
            scale_cmd = 'newstack -float 2 -meansd ' + scale + ' -mode 2 -in ' + inp['invtomo'] + ' -ou ' + inp['scaletomo']
            subprocess.run(scale_cmd, shell=True)
            logf.write("##### I have scaled the density and changed the mode to 32bit floating point with this command: #####'\n'")
            logf.write(scale_cmd + '\n\n')
            
            if decon_para['do_fft'] == 'Yes':
                do_fft(inputfile)
            if decon_para['do_ps'] == 'Yes':
                do_ps(inputfile)                
        else:
            logf.write('I will use the existing inverted and scaled tomogram file: ' + inp['scaletomo'])

        logf.write("##### The header of the scaled tomogram #####'\n'")
        subprocess.run(['header', inp['scaletomo']], stdout=logf)
        logf.write('\n\n')


def do_fft(inputfile):
    with open(inp['logfile'], 'a') as logf:
        logf.write("### I will create a summed FFT with these commands: ###'\n'")
        if not os.path.exists('FFT'):
            os.mkdir('FFT')
        # reorient the files
        if not os.path.isfile('FFT/flipxz_' + inputfile): 
            rot_xz_origtomo_cmd = 'clip flipxz ' + inputfile + ' FFT/flipxz_' + inputfile
            logf.write(rot_xz_origtomo_cmd + '\n\n')
            subprocess.run(rot_xz_origtomo_cmd, shell=True)
        if not os.path.isfile('FFT/flipyz_' + inputfile):
            rot_yz_origtomo_cmd = 'clip flipyz ' + inputfile + ' FFT/flipyz_' + inputfile
            logf.write(rot_yz_origtomo_cmd + '\n\n')
            subprocess.run(rot_yz_origtomo_cmd, shell=True)

        # Calculate a summed FFT
        orig_fft_cmd = 'fftrans -input ' + inputfile + ' -output FFT/2D_FFT_' + inputfile
        avg_orig_fft_cmd = 'clip avg -3d FFT/2D_FFT_' + inputfile + ' avg2D_FFT_' + inputfile

        ### Calculate a summed FFT of the flipxz
        rot_xz_fft_cmd = 'fftrans -input FFT/flipxz_' + inputfile + ' -output FFT/2D_FFT_flipxz_' + inputfile
        avg_xz_fft_cmd = 'clip avg -3d FFT/2D_FFT_flipxz_' + inputfile + ' avg2D_FFT_flipxz_' + inputfile

        ### Calculate a summed FFT of the flipyz
        rot_yz_fft_cmd = 'fftrans -input FFT/flipyz_' + inputfile + ' -output FFT/2D_FFT_flipyz_' + inputfile
        avg_yz_fft_cmd = 'clip avg -3d FFT/2D_FFT_flipyz_' + inputfile + ' avg2D_FFT_flipyz_' + inputfile
        
        cmds = [orig_fft_cmd, avg_orig_fft_cmd, rot_xz_fft_cmd, avg_xz_fft_cmd, rot_yz_fft_cmd, avg_yz_fft_cmd]
        for cmd in cmds:
            logf.write(cmd + '\n')
            subprocess.run(cmd, shell=True)
            logf.write('\n\n')
            
def do_ps(inputfile):
    with open(inp['logfile'], 'a') as logf:
        logf.write("### I will create a summed PowerSpectrum with these commands: ###'\n'")
        if not os.path.exists('FFT'):
            os.mkdir('FFT')
        # reorient the files
        if not os.path.isfile('FFT/flipxz_' + inputfile): 
            rot_xz_origtomo_cmd = 'clip flipxz ' + inputfile + ' FFT/flipxz_' + inputfile
            logf.write(rot_xz_origtomo_cmd + '\n\n')
            subprocess.run(rot_xz_origtomo_cmd, shell=True)
        if not os.path.isfile('FFT/flipyz_' + inputfile):
            rot_yz_origtomo_cmd = 'clip flipyz ' + inputfile + ' FFT/flipyz_' + inputfile
            logf.write(rot_yz_origtomo_cmd + '\n\n')
            subprocess.run(rot_yz_origtomo_cmd, shell=True)

        # Calculate a summed powerspectrum
        orig_ps_cmd = 'clip spectrum ' + inputfile + ' FFT/2D_PS_' + inputfile
        avg_orig_ps_cmd = 'clip avg -3d FFT/2D_PS_' + inputfile + ' avg2D_PS_' + inputfile

        ### Calculate a summed powerspectrum of the flipxz
        rot_xz_ps_cmd = 'clip spectrum FFT/flipxz_' + inputfile + ' FFT/2D_PS_flipxz_' + inputfile
        avg_rot_xz_ps_cmd = 'clip avg -3d FFT/2D_PS_flipxz_' + inputfile + ' avg2D_PS_flipxz_' + inputfile

        ### Calculate a summed powerspectrum of the flipyz
        rot_yz_ps_cmd = 'clip spectrum FFT/flipyz_' + inputfile + ' FFT/2D_PS_flipyz_' + inputfile
        avg_rot_yz_ps_cmd = 'clip avg -3d FFT/2D_PS_flipyz_' + inputfile + ' avg2D_PS_flipyz_' + inputfile

        cmds = [orig_ps_cmd, avg_orig_ps_cmd, rot_xz_ps_cmd, avg_rot_xz_ps_cmd, rot_yz_ps_cmd, avg_rot_yz_ps_cmd]
        for cmd in cmds:
            logf.write(cmd + '\n')
            subprocess.run(cmd, shell=True)
            logf.write('\n\n')        
            
def decon():
    with open(inp['logfile'], 'a') as logf:
        logf.write('### I will run deconvolutionruns with' + '\n')
        logf.write('smothing factors of ' + str(decon_para['smooth']) + ' and ' + str(
            decon_para['ni']) + ' iterations ###' + '\n')
        for it in decon_para['ni']:
            for smooth in decon_para['smooth']:
                for lamsmooth in decon_para['lamsmooth']:
                    logf.write("##### NOW COMES THE REAL DECONVOLUTION STEP: #####'\n'")
                    decon_out = 'decon-s' + str(smooth) + '_nc' + str(it) + '_lams' + str(lamsmooth) + '-' + inp['scaletomo']
                    plotfile = 'plotfile-s' + str(smooth) + '_nc' + str(it) + '_lams' + str(lamsmooth) + '-' + inp['tomoname'] + '_inv_scaled'
                    if not os.path.isfile(decon_out):
                        try:
                            decon_cmd = 'core2_decon ' + inp['scaletomo'] + ' ' + decon_out + ' ' + inp['OTF_file'] + ' -alpha=100 -lamratio=0:1 -lamf=' + str(smooth) + ' -lampc=0 -lampos=1 -lamsmooth=' + str(lamsmooth) + ' -laml2=0 -cuth=0.001 -ncycl=' + str(it) + ' -omega=0.8 -sub=1:1:1:1:1 -tol=0.0001 -oplotfile=' + plotfile + ' -linesearch="2014" -regtype="ma" -logging=progress' #-np=' + decon_para['np'] + '
                            logf.write(decon_cmd + '\n\n')
                            subprocess.run(decon_cmd, shell=True, stdout=logf)
                            logf.write('\n\n')
                            deconheader = subprocess.run(['header', decon_out], capture_output=True, text=True)  # read the header into stdout
                            logf.write(deconheader.stdout)  # write the header into the logfile
                            histo(decon_out)
                        except:
                            print('############ Decon failed! ############')
                            print('Maybe because the 3D volume was too big.')
                            print('Try to bin your data!')
                    else:
                        logf.write('There is already a file called ' + decon_out)
                    if decon_para['do_fft'] == 'Yes':
                        do_fft(decon_out)
                    if decon_para['do_ps'] == 'Yes':
                        do_ps(decon_out)

def histo(in_file):
    with open(inp['logfile'], 'a') as logf:
        logf.write('### I will scale the histogram of ' + in_file + 'with this command: \n\n')  
        percentile_cmd = 'mrcpercentiles ' + in_file + ' -cut=0.001 -cut=0.999'
        percentile = subprocess.run(percentile_cmd, shell=True, capture_output=True, text=True)
        logf.write(percentile_cmd)
        logf.write(percentile.stdout)
        perc_out = subprocess.run(['grep', test_dict['xy']], capture_output=True, text=True, input=percentile.stdout)
        low = perc_out.stdout.split()[-2]
        high = perc_out.stdout.split()[-1]
        truncate_cmd = 'clip truncate -l ' + str(low) + ' -h ' + str(high) + ' ' + in_file + ' ' + in_file
        logf.write(truncate_cmd)
        subprocess.run(truncate_cmd, shell=True)
        if os.path.exists(in_file + '~'):
            os.remove(in_file + '~')

        
def clean_up(): 
    #And delete intermediate files
    trees = ['FFT', 'rotvol', 'rotvol_b']
    removers = [inp['invtomo'], inp['scaletomo']]
    for tree in trees:
        if os.path.exists(tree):
            shutil.rmtree(tree)
    for remover in removers:
        if os.path.exists(remover):
            print (remover)
            os.remove(remover)
    #Now I clear the dictionaries for the next run
    decon_para.clear()
    inp.clear()
    test_dict.clear()



   
#######################
### Main Parameters ###
#######################

#homedir = '/home/labs/peterk/peterkir/'
parentdir = os.getcwd()  # directory were the script and original tomograms are located
print(parentdir)
environment()
scale = '0,1000'
#scale = '30000,2000'

for mrc_file in os.listdir('.'):
    if fnmatch.fnmatch(mrc_file, '*_rec.mrc'):
        start_time = perf_counter()
        input_mrc = mrc_file
        inp = input_files(input_mrc, parentdir)
        write_log_file(inp['logfile'])
        test_dict = extract_data(inp['newfile'])
        decon_para = de_para()
        prepare_input_files(input_mrc)
        fan_prep()
        ftransform()
        mid_time = perf_counter()
        decon()
        stop_time = perf_counter()
        fan_time = mid_time-start_time
        used_time = stop_time-start_time
        with open(inp['logfile'], 'a') as logf: 
            logf.write(f'FanOut finished in {round(fan_time,2)} seconds, {round(fan_time/60,2)} minutes, or {round(fan_time/3600,2)} hours'+ '\n')
            logf.write(f'Whole run finished in {round(used_time,2)} seconds, {round(used_time/60,2)} minutes, or {round(used_time/3600,2)} hours')
        clean_up()

