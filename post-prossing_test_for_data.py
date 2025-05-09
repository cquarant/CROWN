import argparse
import sys
import yaml
import os
import ROOT as R
import time
import array
import socket

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

def post_proc_varial(input_path, output_path, samples_list):
    list_all = getall_list(input_path) 
    for f in list_all:
        print("found")
        try:
            df_mc = R.RDataFrame('ntuple', f)
            print("define df")
            df_mc = df_mc.Filter('(HLT_AK8PFHT800_TrimMass50 || HLT_AK8PFJet400_TrimMass30 || HLT_AK8PFJet500 ||HLT_PFJet500||HLT_PFHT1050||HLT_PFHT500_PFMET100_PFMHT100_IDTight||HLT_PFHT700_PFMET85_PFMHT85_IDTight||HLT_PFHT800_PFMET75_PFMHT75_IDTight)')
            print('hlt finish')
            #col_names = df_mc.GetColumnNames()
        except:
            print("No Column get in this file {}".format(f))
            continue
        #if 'Xsec'  not in col_names:
        df_mc = df_mc.Define('Xsec', '1.0').Define('genEventSumW', '1.0') 
        print("xs and gw added")
        df_mc = df_mc.Define('Train_weight','1.0')
        print("tw added")
        df_mc = df_mc.Define('puweight','1.0')
        print("pw added")
        df_mc.Snapshot('ntuple',  f.replace(input_path, output_path))
        print("this channnel finished")

if __name__ == '__main__':       

    input_path = "/eos/user/j/jinwa/CROWN_output_first_dataMC_backup"
    output_path = "/eos/user/j/jinwa/CROWN_post_output"

    with open("/afs/cern.ch/user/j/jinwa/KingMaker/sample_database/datasets.yaml" , "r") as file:
        samples_list =  yaml.safe_load(file)
        #print(samples_list)
        post_proc_varial(input_path, output_path, samples_list)
