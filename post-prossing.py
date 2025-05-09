import argparse
import sys
import yaml
import os
import ROOT as R
import time
import array
import socket

R.gROOT.SetBatch(True)
R.gInterpreter.Declare('#include "/data/bond/botaoguo/CROWN/build_run3/bin/applyfactor.h"')

# need CROWN env
lumi_dict = {
    # https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmVRun3Analysis
    # https://twiki.cern.ch/twiki/bin/view/CMS/BrilcalcQuickStart
    # Brilcalc # 63.145882859
    # 34.748619391999995
    "2022preEE" : 8.077009685e3, # pb^-1
    "2022postEE": 26.671609707e3,# pb^-1 
    # 28.397263467000002
    "2023preBPix": 18.704133414e3, # B+C
    "2023postBPix": 9.693130053e3,
}

# 2022
fake_factor_dict = {
    "2022preEE": 
        {
            "e2m": 0.022747415392360283,
            "m2m": 0.024111873284901615,
        },
    "2022postEE": 
        {
            "e2m": 0.022747415392360283,
            "m2m": 0.024111873284901615,
        },
    "2023preBPix":
        {
            "e2m": 0.024995812473941486,
            "m2m": 0.024334559187637824,
        },
    "2023postBPix":
        {
            "e2m": 0.024995812473941486,
            "m2m": 0.024334559187637824,
        },    
}

genWeight_dict = {
    # preEE  abs(genWeight) < 3
    'WminusH_Hto2Mu_WtoAll_M-125_TuneCP5_13p6TeV_powheg-pythia8_Run3Summer22NanoAODv12-130X' : 77443.426,
    'WplusH_Hto2Mu_WtoAll_M-125_TuneCP5_13p6TeV_powheg-pythia8_Run3Summer22NanoAODv12-130X'  : 112512.06,
    'ZH_Hto2Mu_ZtoAll_M-125_TuneCP5_13p6TeV_powheg-pythia8_Run3Summer22NanoAODv12-130X'      : 108873.63,
    # postEE  abs(genWeight) < 3
    'WminusH_Hto2Mu_WtoAll_M-125_TuneCP5_13p6TeV_powheg-pythia8_Run3Summer22EENanoAODv12-130X' : 271691.31,
    'WplusH_Hto2Mu_WtoAll_M-125_TuneCP5_13p6TeV_powheg-pythia8_Run3Summer22EENanoAODv12-130X'  : 397696.12,
    'ZH_Hto2Mu_ZtoAll_M-125_TuneCP5_13p6TeV_powheg-pythia8_Run3Summer22EENanoAODv12-130X'      : 380191.51,
    # preBPix  abs(genWeight) < 3
    'WminusH_Hto2Mu_WtoAll_M-125_TuneCP5_13p6TeV_powheg-minlo-pythia8_Run3Summer23NanoAODv12-130X' : 234383.14,
    'WplusH_Hto2Mu_WtoAll_M-125_TuneCP5_13p6TeV_powheg-minlo-pythia8_Run3Summer23NanoAODv12-130X'  : 340136.73,
    'ZH_Hto2Mu_ZtoAll_M-125_TuneCP5_13p6TeV_powheg-minlo-pythia8_Run3Summer23NanoAODv12-130X'      : 324608.83,
    # postBPix abs(genWeight) < 3
    'WminusH_Hto2Mu_WtoAll_M-125_TuneCP5_13p6TeV_powheg-minlo-pythia8_Run3Summer23BPixNanoAODv12-130X' : 74386.730,
    'WplusH_Hto2Mu_WtoAll_M-125_TuneCP5_13p6TeV_powheg-minlo-pythia8_Run3Summer23BPixNanoAODv12-130X'  : 113017.93,
    'ZH_Hto2Mu_ZtoAll_M-125_TuneCP5_13p6TeV_powheg-minlo-pythia8_Run3Summer23BPixNanoAODv12-130X'      : 164743.67,
}

col_dict = { 
    'Train_weight': '1.0',
    'puweight': '1.0',
    'puweight__PileUpDown': '1.0',
    'puweight__PileUpUp' : '1.0',
    'BosonDecayMode': '-10',
    
    'btag_weight': '1.0f',
    'genWeight': '1.0f', 
    'Xsec' : '1.0f', 
    'genEventSumW': '1.0f', 
    'genmet_pt': '-10.f',
    'genmet_phi': '-10.f',
    
}
#
process_name = [
    'WplusH_Hto2Mu','WminusH_Hto2Mu','ZH_Hto2Mu', 'TTH_Hto2Mu', "GluGluHto2Mu", "VBFHto2Mu", # vhmm
    'WZto3LNu','WZto2L2Q','WWto2L2Nu','ZZto2L2Q','ZZto2L2Nu','ZZto4L', # diboson
    'TTto2L2Nu','TWminusto2L2Nu','TbarWplusto2L2Nu', # top
    # dyjets and data only 1 process, no need to add flag
    'WWW','WWZ','WZZ','ZZZ', # triboson
]
# add is_{process_name} to col_dict
for id in process_name:
    col_dict.update( {'is_'+id: '0.0f'} ) 

def ListToStdVector(l, elem_type='string'):
    v = R.std.vector(elem_type)()
    for x in l:
        if elem_type in ['Int_t', 'UInt_t']:
            x = int(x)
        v.push_back(x)
    return v


def getall_list(input_path):
    '''
    get a list of all files in input_path
    '''
    inputfile = []
    for fname in os.listdir(input_path):
        if '.root' not in fname or 'FakeFactor' in fname or 'output' in fname or 'input' in fname or 'Fakes' in fname or 'FF' in fname:
            continue
        _file = R.TFile.Open(input_path + '/' + fname)
        if _file.GetListOfKeys().Contains("ntuple"):
            inputfile.append(input_path + '/' + fname)
    # print(inputfile)
    # exit(0)
    return inputfile

def Add_new_column(input_path, output_path, samples_list, new_column):
    '''
    add a new column, new_column should be a dictionary with key name of 
    new column with value of column value (check what type you need)
    new_column = {column_name: value}
    '''
    f_list = getall_list(input_path)
    for f in f_list:
        # add weight==1 only for data
        if 'SingleMuon' not in f and 'Muon' not in f:
            continue
        # loop over each file 
        # if update regionc weight
        if Updated_regionc:
            if "e2m_dyfakeinge_regionc" in f:
                new_column.update({'Train_weight': fake_factor_dict[year]["e2m"]})
            elif "m2m_dyfakeingmu_regionc" in f:
                new_column.update({'Train_weight': fake_factor_dict[year]["m2m"]})
        for n in samples_list:
            if n in f:
                print(f)
                try:
                    df_mc = R.RDataFrame('ntuple', f)
                except:
                    print("no tree named ntuple in this file")
                    continue
                col_names = df_mc.GetColumnNames()
                if len(col_names) <=0:
                    print("{} has no branch".format(f))
                    continue
                # print(col_names)
                # exit(0)
                edited = False
                for c in new_column:    
                    if c not in col_names:    
                        edited = True
                        df_mc = df_mc.Define(c,str(new_column[c]))
                if "e2m" in f  or "m2m" in f:
                    df_mc = df_mc.Define('lep_MET_dphi', "(extra_lep_phi - met_phi > M_PI) ? (extra_lep_phi - met_phi - 2.0*M_PI) : ( (extra_lep_phi - met_phi <= - M_PI) ? (extra_lep_phi - met_phi + 2.0*M_PI) : (extra_lep_phi - met_phi) )")
                    if "2022" in year:
                        df_mc = df_mc.Define("wz_zz_scale","applyWZscale2022(is_diboson, is_WWto2L2Nu, nelectrons, nbaseelectrons, nmuons, nbasemuons)")
                    elif "2023" in year:
                        df_mc = df_mc.Define("wz_zz_scale","applyWZscale2023(is_diboson, is_WWto2L2Nu, nelectrons, nbaseelectrons, nmuons, nbasemuons)")
                # commom sfs
                df_mc = df_mc.Define("id_wgt_mu_1","1.0")
                df_mc = df_mc.Define("iso_wgt_mu_1","1.0")
                df_mc = df_mc.Define("id_wgt_mu_1_below15","1.0")
                df_mc = df_mc.Define("iso_wgt_mu_1_below15","1.0")
                df_mc = df_mc.Define("id_wgt_mu_2","1.0")
                df_mc = df_mc.Define("iso_wgt_mu_2","1.0")
                df_mc = df_mc.Define("id_wgt_mu_2_below15","1.0")
                df_mc = df_mc.Define("iso_wgt_mu_2_below15","1.0")
                
                df_mc = df_mc.Define("id_wgt_mu_1__MuonIDIsoUp","1.0")
                df_mc = df_mc.Define("iso_wgt_mu_1__MuonIDIsoUp","1.0")
                df_mc = df_mc.Define("id_wgt_mu_1_below15__MuonIDIsoUp","1.0")
                df_mc = df_mc.Define("iso_wgt_mu_1_below15__MuonIDIsoUp","1.0")
                df_mc = df_mc.Define("id_wgt_mu_2__MuonIDIsoUp","1.0")
                df_mc = df_mc.Define("iso_wgt_mu_2__MuonIDIsoUp","1.0")
                df_mc = df_mc.Define("id_wgt_mu_2_below15__MuonIDIsoUp","1.0")
                df_mc = df_mc.Define("iso_wgt_mu_2_below15__MuonIDIsoUp","1.0")

                df_mc = df_mc.Define("id_wgt_mu_1__MuonIDIsoDown","1.0")
                df_mc = df_mc.Define("iso_wgt_mu_1__MuonIDIsoDown","1.0")
                df_mc = df_mc.Define("id_wgt_mu_1_below15__MuonIDIsoDown","1.0")
                df_mc = df_mc.Define("iso_wgt_mu_1_below15__MuonIDIsoDown","1.0")
                df_mc = df_mc.Define("id_wgt_mu_2__MuonIDIsoDown","1.0")
                df_mc = df_mc.Define("iso_wgt_mu_2__MuonIDIsoDown","1.0")
                df_mc = df_mc.Define("id_wgt_mu_2_below15__MuonIDIsoDown","1.0")
                df_mc = df_mc.Define("iso_wgt_mu_2_below15__MuonIDIsoDown","1.0")

                # 3l, 4l
                if "e2m" in f  or "m2m" in f or "eemm" in f or "mmmm" in f:
                    df_mc = df_mc.Define("id_wgt_mu_3","1.0")
                    df_mc = df_mc.Define("iso_wgt_mu_3","1.0")
                    df_mc = df_mc.Define("id_wgt_mu_3_below15","1.0")
                    df_mc = df_mc.Define("iso_wgt_mu_3_below15","1.0")

                    df_mc = df_mc.Define("id_wgt_ele_loose_1","1.0")
                    df_mc = df_mc.Define("id_wgt_ele_wp90Iso_1","1.0")

                    df_mc = df_mc.Define("id_wgt_mu_3__MuonIDIsoUp","1.0")
                    df_mc = df_mc.Define("iso_wgt_mu_3__MuonIDIsoUp","1.0")
                    df_mc = df_mc.Define("id_wgt_mu_3_below15__MuonIDIsoUp","1.0")
                    df_mc = df_mc.Define("iso_wgt_mu_3_below15__MuonIDIsoUp","1.0")

                    df_mc = df_mc.Define("id_wgt_mu_3__MuonIDIsoDown","1.0")
                    df_mc = df_mc.Define("iso_wgt_mu_3__MuonIDIsoDown","1.0")
                    df_mc = df_mc.Define("id_wgt_mu_3_below15__MuonIDIsoDown","1.0")
                    df_mc = df_mc.Define("iso_wgt_mu_3_below15__MuonIDIsoDown","1.0")

                    df_mc = df_mc.Define("id_wgt_ele_loose_1__EleIDUp","1.0")
                    df_mc = df_mc.Define("id_wgt_ele_wp90Iso_1__EleIDUp","1.0")
                    
                    df_mc = df_mc.Define("id_wgt_ele_loose_1__EleIDDown","1.0")
                    df_mc = df_mc.Define("id_wgt_ele_wp90Iso_1__EleIDDown","1.0")
                    
                    df_mc = df_mc.Define("reco_wgt_ele_1","1.0")
                    df_mc = df_mc.Define("reco_wgt_ele_1__EleRecoDown","1.0")
                    df_mc = df_mc.Define("reco_wgt_ele_1__EleRecoUp","1.0")

                # 4l
                if "eemm" in f or "mmmm" in f:
                    df_mc = df_mc.Define("id_wgt_mu_4","1.0")
                    df_mc = df_mc.Define("iso_wgt_mu_4","1.0")
                    df_mc = df_mc.Define("id_wgt_mu_4_below15","1.0")
                    df_mc = df_mc.Define("iso_wgt_mu_4_below15","1.0")

                    df_mc = df_mc.Define("id_wgt_ele_loose_2","1.0")
                    df_mc = df_mc.Define("id_wgt_ele_wp90Iso_2","1.0")

                    df_mc = df_mc.Define("id_wgt_mu_4__MuonIDIsoUp","1.0")
                    df_mc = df_mc.Define("iso_wgt_mu_4__MuonIDIsoUp","1.0")
                    df_mc = df_mc.Define("id_wgt_mu_4_below15__MuonIDIsoUp","1.0")
                    df_mc = df_mc.Define("iso_wgt_mu_4_below15__MuonIDIsoUp","1.0")

                    df_mc = df_mc.Define("id_wgt_mu_4__MuonIDIsoDown","1.0")
                    df_mc = df_mc.Define("iso_wgt_mu_4__MuonIDIsoDown","1.0")
                    df_mc = df_mc.Define("id_wgt_mu_4_below15__MuonIDIsoDown","1.0")
                    df_mc = df_mc.Define("iso_wgt_mu_4_below15__MuonIDIsoDown","1.0")

                    df_mc = df_mc.Define("id_wgt_ele_loose_2__EleIDUp","1.0")
                    df_mc = df_mc.Define("id_wgt_ele_wp90Iso_2__EleIDUp","1.0")

                    df_mc = df_mc.Define("id_wgt_ele_loose_2__EleIDDown","1.0")
                    df_mc = df_mc.Define("id_wgt_ele_wp90Iso_2__EleIDDown","1.0")
                    
                    df_mc = df_mc.Define("reco_wgt_ele_2","1.0")
                    df_mc = df_mc.Define("reco_wgt_ele_2__EleRecoDown","1.0")
                    df_mc = df_mc.Define("reco_wgt_ele_2__EleRecoUp","1.0")

                if "nnmm" in f:
                    # add gen mu 
                    df_mc = df_mc.Define("genmu1_fromH_pt","-1.0f")
                    df_mc = df_mc.Define("genmu1_fromH_eta","-1.0f")
                    df_mc = df_mc.Define("genmu1_fromH_phi","-1.0f")
                    df_mc = df_mc.Define("genmu1_fromH_mass","-1.0f")
                    df_mc = df_mc.Define("genmu2_fromH_pt","-1.0f")
                    df_mc = df_mc.Define("genmu2_fromH_eta","-1.0f")
                    df_mc = df_mc.Define("genmu2_fromH_phi","-1.0f")
                    df_mc = df_mc.Define("genmu2_fromH_mass","-1.0f")
                # apply DY scale in fjmm
                if "fjmm" in f:
                    df_mc = df_mc.Define("dy_scale", "1.0f")
                if "2016" in year:
                    df_mc = df_mc.Define("trg_single_mu27", "true")
                    df_mc = df_mc.Define("is_2016", "1.0f")
                else:
                    df_mc = df_mc.Define("is_2016", "0.0f")
                if "2017" in year:
                    df_mc = df_mc.Define("is_2017", "1.0f")
                else:
                    df_mc = df_mc.Define("is_2017", "0.0f")
                if "2018" in year:
                    df_mc = df_mc.Define("is_2018", "1.0f")
                else:
                    df_mc = df_mc.Define("is_2018", "0.0f")                
                if "2022" in year:
                    df_mc = df_mc.Define("is_2022", "1.0f")
                else:
                    df_mc = df_mc.Define("is_2022", "0.0f")
                if "2023" in year:
                    df_mc = df_mc.Define("is_2023", "1.0f")
                else:
                    df_mc = df_mc.Define("is_2023", "0.0f")
                
                col_names_v2 = df_mc.GetColumnNames()
                output_list = []
                for i in range(len(col_names_v2)):
                    if "_corrected" in col_names_v2[i]:
                        continue
                    output_list.append(col_names_v2[i])
                # for i in new_column.keys():
                #     output_list.append(i)

                if edited:
                    df_mc.Snapshot('ntuple',  f.replace(input_path, output_path), ListToStdVector(output_list))
                else:
                    print( 'column {} already exists, skip'.format(new_column) )

def post_proc_varial(input_path, output_path, samples_list):
    list_all = getall_list(input_path) 
    for f in list_all:
        if 'SingleMuon' in f or 'Muon' in f:
            continue
        # remove the sig in CR
        if "cr" in channel:
            if "Hto2Mu" in f:
                continue
        for n in samples_list:
            if n in f:
                # print(n)
                print(f)
                try:
                    df_mc = R.RDataFrame('ntuple', f)
                    col_names = df_mc.GetColumnNames()
                except:
                    print("No Column get in this file {}".format(f))
                    continue
                for id in process_name:
                    if id in f:
                        df_mc = df_mc.Define('is_'+id, str('1.0f'))
                    else:
                        df_mc = df_mc.Define('is_'+id, str('0.0f'))
                ###############################################
                # create an empty flag for mc in e2m region
                # if 'e2m_dyfakeinge_regionc' in f:
                #     df_mc = df_mc.Define('Flag_dimuon_Zmass_veto', str('1.0f'))
                ###############################################
                if 'Xsec'  not in col_names:
                    # cut the abnormal region in 2022postEE
                    # please cut event that has the abnormal genWeight (currently only for vhmm signal?)
                    if 'WplusH_Hto2Mu' in f or 'WminusH_Hto2Mu' in f or 'ZH_Hto2Mu' in f:
                        df_mc = df_mc.Filter('abs(genWeight)<3')
                        # need to change the genWeightSum of signal
                        gensumw = genWeight_dict[n]
                        print("abnormal genweight event processing at {} !!!!!!!!!!".format(n))
                    else:
                        # if merge e2m and m2m, plz use  gensumw /= 2
                        # if only single channel, e2m or m2m, use gensumw
                        gensumw = R.RDataFrame('conditions', f).Sum('genEventSumw').GetValue()
                        # print("gensumw : {}".format(gensumw))
                    if "Embedding" in f:
                        df_mc = df_mc.Define('Xsec', str(samples_list[n]['xsec']) +'f').Define('genEventSumW', '1.0f') # set gensumw for embedding as 1, since we use genWeight directly
                    else:
                        df_mc = df_mc.Define('Xsec', str(samples_list[n]['xsec']) +'f').Define('genEventSumW', str(gensumw) + 'f')
                    # df_mc = df_mc.Define('Train_weight','Xsec* 139e3 * genWeight / genEventSumW') ## define weight calculated
                    lumi = lumi_dict[year]
                    if Updated_regionc:
                        if "e2m_dyfakeinge_regionc" in f:
                            df_mc = df_mc.Define('Train_weight','( -{0} * Xsec* {1} * genWeight ) / genEventSumW'.format(fake_factor_dict[year]["e2m"], lumi_dict[year]))
                        elif "m2m_dyfakeingmu_regionc" in f:
                            df_mc = df_mc.Define('Train_weight','( -{0} * Xsec* {1} * genWeight ) / genEventSumW'.format(fake_factor_dict[year]["m2m"], lumi_dict[year]))
                    else:
                        df_mc = df_mc.Define('Train_weight','( 1.0 * Xsec* {} * genWeight ) / genEventSumW'.format(lumi_dict[year]))
                    if "e2m" in f or "m2m" in f:
                        df_mc = df_mc.Define('lep_MET_dphi', "(extra_lep_phi - met_phi > M_PI) ? (extra_lep_phi - met_phi - 2.0*M_PI) : ( (extra_lep_phi - met_phi <= - M_PI) ? (extra_lep_phi - met_phi + 2.0*M_PI) : (extra_lep_phi - met_phi) )")
                        if "202" in year:
                            if "e2m" in f:
                                df_mc = df_mc.Define("id_wgt_mu_3","1.0")
                                df_mc = df_mc.Define("iso_wgt_mu_3","1.0")
                                df_mc = df_mc.Define("id_wgt_mu_3_below15","1.0")
                                df_mc = df_mc.Define("iso_wgt_mu_3_below15","1.0")
                                
                                df_mc = df_mc.Define("id_wgt_mu_3__MuonIDIsoUp","1.0")
                                df_mc = df_mc.Define("iso_wgt_mu_3__MuonIDIsoUp","1.0")
                                df_mc = df_mc.Define("id_wgt_mu_3_below15__MuonIDIsoUp","1.0")
                                df_mc = df_mc.Define("iso_wgt_mu_3_below15__MuonIDIsoUp","1.0")
                                
                                df_mc = df_mc.Define("id_wgt_mu_3__MuonIDIsoDown","1.0")
                                df_mc = df_mc.Define("iso_wgt_mu_3__MuonIDIsoDown","1.0")
                                df_mc = df_mc.Define("id_wgt_mu_3_below15__MuonIDIsoDown","1.0")
                                df_mc = df_mc.Define("iso_wgt_mu_3_below15__MuonIDIsoDown","1.0")

                            elif "m2m" in f:
                                df_mc = df_mc.Define("id_wgt_ele_loose_1","1.0")
                                df_mc = df_mc.Define("id_wgt_ele_wp90Iso_1","1.0")
                                df_mc = df_mc.Define("id_wgt_ele_loose_1__EleIDUp","1.0")
                                df_mc = df_mc.Define("id_wgt_ele_wp90Iso_1__EleIDUp","1.0")
                                df_mc = df_mc.Define("id_wgt_ele_loose_1__EleIDDown","1.0")
                                df_mc = df_mc.Define("id_wgt_ele_wp90Iso_1__EleIDDown","1.0")

                                df_mc = df_mc.Define("reco_wgt_ele_1","1.0")
                                df_mc = df_mc.Define("reco_wgt_ele_1__EleRecoDown","1.0")
                                df_mc = df_mc.Define("reco_wgt_ele_1__EleRecoUp","1.0")

                            if "2022" in year:
                                df_mc = df_mc.Define("wz_zz_scale","applyWZscale2022(is_diboson, is_WWto2L2Nu, nelectrons, nbaseelectrons, nmuons, nbasemuons)")
                            elif "2023" in year:
                                df_mc = df_mc.Define("wz_zz_scale","applyWZscale2023(is_diboson, is_WWto2L2Nu, nelectrons, nbaseelectrons, nmuons, nbasemuons)")
                    if "eemm" in f:
                        df_mc = df_mc.Define("id_wgt_mu_3","1.0")
                        df_mc = df_mc.Define("iso_wgt_mu_3","1.0")
                        df_mc = df_mc.Define("id_wgt_mu_3_below15","1.0")
                        df_mc = df_mc.Define("iso_wgt_mu_3_below15","1.0")

                        df_mc = df_mc.Define("id_wgt_mu_4","1.0")
                        df_mc = df_mc.Define("iso_wgt_mu_4","1.0")
                        df_mc = df_mc.Define("id_wgt_mu_4_below15","1.0")
                        df_mc = df_mc.Define("iso_wgt_mu_4_below15","1.0")

                        df_mc = df_mc.Define("id_wgt_mu_3__MuonIDIsoUp","1.0")
                        df_mc = df_mc.Define("iso_wgt_mu_3__MuonIDIsoUp","1.0")
                        df_mc = df_mc.Define("id_wgt_mu_3_below15__MuonIDIsoUp","1.0")
                        df_mc = df_mc.Define("iso_wgt_mu_3_below15__MuonIDIsoUp","1.0")
                        df_mc = df_mc.Define("id_wgt_mu_3__MuonIDIsoDown","1.0")
                        df_mc = df_mc.Define("iso_wgt_mu_3__MuonIDIsoDown","1.0")
                        df_mc = df_mc.Define("id_wgt_mu_3_below15__MuonIDIsoDown","1.0")
                        df_mc = df_mc.Define("iso_wgt_mu_3_below15__MuonIDIsoDown","1.0")

                        df_mc = df_mc.Define("id_wgt_mu_4__MuonIDIsoUp","1.0")
                        df_mc = df_mc.Define("iso_wgt_mu_4__MuonIDIsoUp","1.0")
                        df_mc = df_mc.Define("id_wgt_mu_4_below15__MuonIDIsoUp","1.0")
                        df_mc = df_mc.Define("iso_wgt_mu_4_below15__MuonIDIsoUp","1.0")
                        df_mc = df_mc.Define("id_wgt_mu_4__MuonIDIsoDown","1.0")
                        df_mc = df_mc.Define("iso_wgt_mu_4__MuonIDIsoDown","1.0")
                        df_mc = df_mc.Define("id_wgt_mu_4_below15__MuonIDIsoDown","1.0")
                        df_mc = df_mc.Define("iso_wgt_mu_4_below15__MuonIDIsoDown","1.0")

                    elif "mmmm" in f:
                        df_mc = df_mc.Define("id_wgt_ele_loose_1","1.0")
                        df_mc = df_mc.Define("id_wgt_ele_wp90Iso_1","1.0")
                        df_mc = df_mc.Define("id_wgt_ele_loose_2","1.0")
                        df_mc = df_mc.Define("id_wgt_ele_wp90Iso_2","1.0")

                        df_mc = df_mc.Define("id_wgt_ele_loose_1__EleIDUp","1.0")
                        df_mc = df_mc.Define("id_wgt_ele_wp90Iso_1__EleIDUp","1.0")
                        df_mc = df_mc.Define("id_wgt_ele_loose_1__EleIDDown","1.0")
                        df_mc = df_mc.Define("id_wgt_ele_wp90Iso_1__EleIDDown","1.0")
                        df_mc = df_mc.Define("id_wgt_ele_loose_2__EleIDUp","1.0")
                        df_mc = df_mc.Define("id_wgt_ele_wp90Iso_2__EleIDUp","1.0")
                        df_mc = df_mc.Define("id_wgt_ele_loose_2__EleIDDown","1.0")
                        df_mc = df_mc.Define("id_wgt_ele_wp90Iso_2__EleIDDown","1.0")

                        df_mc = df_mc.Define("reco_wgt_ele_1","1.0")
                        df_mc = df_mc.Define("reco_wgt_ele_1__EleRecoDown","1.0")
                        df_mc = df_mc.Define("reco_wgt_ele_1__EleRecoUp","1.0")
                        df_mc = df_mc.Define("reco_wgt_ele_2","1.0")
                        df_mc = df_mc.Define("reco_wgt_ele_2__EleRecoDown","1.0")
                        df_mc = df_mc.Define("reco_wgt_ele_2__EleRecoUp","1.0")

                    # apply dy scale in fjmm
                    if "fjmm" in f and "2022" in year:
                        df_mc = df_mc.Define("dy_scale", "applyDYscale2022(is_dyjets)")
                    elif "fjmm" in f and "2023" in year:
                        df_mc = df_mc.Define("dy_scale", "applyDYscale2023(is_dyjets)")
                    ### add year in var
                    if "2016" in year:
                        df_mc = df_mc.Define("trg_single_mu27", "true")
                        df_mc = df_mc.Define("is_2016", "1.0f")
                    else:
                        df_mc = df_mc.Define("is_2016", "0.0f")
                    if "2017" in year:
                        df_mc = df_mc.Define("is_2017", "1.0f")
                    else:
                        df_mc = df_mc.Define("is_2017", "0.0f")
                    if "2018" in year:
                        df_mc = df_mc.Define("is_2018", "1.0f")
                    else:
                        df_mc = df_mc.Define("is_2018", "0.0f")
                    if "2022" in year:
                        df_mc = df_mc.Define("is_2022", "1.0f")
                    else:
                        df_mc = df_mc.Define("is_2022", "0.0f")
                    if "2023" in year:
                        df_mc = df_mc.Define("is_2023", "1.0f")
                    else:
                        df_mc = df_mc.Define("is_2023", "0.0f")
                    df_mc.Snapshot('ntuple',  f.replace(input_path, output_path)) # update the file finished
    Add_new_column(input_path, output_path, samples_list, col_dict)

if __name__ == '__main__':       
    
    parser = argparse.ArgumentParser(description='post-process for hmm')
    parser.add_argument('--era', required=True, type=str, help="2022preEE, 2022postEE, 2023preBPix, 2023postBPix")
    parser.add_argument('--ch', required=True, type=str, help="2l,3l,4l,cr_bd,cr_c,cr_fjmm")
    parser.add_argument('--updated', required=False, type=bool, default=False, help="update the factor on CR")
    parser.add_argument('--updatedC', required=False, type=bool, default=False, help="update the weight at CR_regionC")
    args = parser.parse_args()

    Updated_regionc = args.updatedC
    year = args.era
    channel = args.ch
    
    current_path = os.getcwd()
    input_path = current_path + '/' + 'input_test_{}'.format(channel)
    if args.updated == False:
        output_path = current_path + '/' + 'output_test_{0}_{1}'.format(channel,year)
    elif args.updated == True:
        output_path = current_path + '/' + 'output_test_{0}_{1}_updated'.format(channel,year)
    os.system("mkdir -p {}".format(output_path))
    
    with open("/data/bond/botaoguo/CROWN/sample_database/datasets{}.yaml".format(year) , "r") as file:
        samples_list =  yaml.safe_load(file)
        # print(samples_list)
        post_proc_varial(input_path, output_path, samples_list)
    print("Target path: {}".format(output_path))
    