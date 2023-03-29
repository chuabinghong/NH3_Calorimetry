import os

## FILE DIRECTORY SETTINGS

wd = r"../"
os.chdir(wd)

data_dir = r"i_data_processed/"
hf_dir = r"o_heatFlow/"

## SCRIPT ARGUMENTS
peak_arg = 0 #1 if we are getting peaks
range_arg = 1 #1 if using -196C to 46C
ramp_rate = 0.1

# adjust experimental data according to range_arg
if range_arg == 0:
    exp_dir = r"i_data/0.10Kmin/"
    mass_ls = [2.5108,4.6293,4.5858,4.5202,3.7778]
    wt_ls = [0, 5.2, 8.4, 10.0, 26.93]  # based off liquidus alignment offset
    t_end = 12981

elif range_arg == 1:
    exp_dir = r"i_data/46C/"
    mass_ls = [4.5386,4.1943,3.8153,3.7107]
    wt_ls = [5.2, 8.2, 14.4, 20.1]  # based off liquidus alignment
    t_end = 13066


os.chdir(r'C:\1_data\OneDrive - Nanyang Technological University\OFYP\CalorimetryAnalysis\z_scripts')
exec(open('1a_hf_prep.py').read())
os.chdir(r'C:\1_data\OneDrive - Nanyang Technological University\OFYP\CalorimetryAnalysis\z_scripts')
exec(open('2a_cp_calc.py').read())
os.chdir(r'C:\1_data\OneDrive - Nanyang Technological University\OFYP\CalorimetryAnalysis\z_scripts')
exec(open('2b_cp_mean.py').read())
os.chdir(r'C:\1_data\OneDrive - Nanyang Technological University\OFYP\CalorimetryAnalysis\z_scripts')
# exec(open('3a_cp_fits.py').read())


