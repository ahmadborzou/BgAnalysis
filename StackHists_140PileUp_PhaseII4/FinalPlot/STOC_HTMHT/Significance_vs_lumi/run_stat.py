#!/usr/bin/env python


signal=4381.63    ## yield signal
background=17752  ## yield background
unc_const=0.0359745*background  ##  rel uncert constant over lumi / abs uncert scale wiht L 
unc_lumi=0.0203829*background  ## uncert scale with 1/sqrt(L) / abs uncert scale with sqrt(L) 
lumi=3000  ## corresponding lumi 


#lumis=[100,3000]
lumis=[100,200,300,500,700,1000,1500,2000,2500,3000]



import sys,ConfigParser,os,string,commands,time,shutil,subprocess,math
from ROOT import TMath,RooStats


class SingleBin:
  def __init__(self,signal,background,constuncert,lumiuncert,lumi):
    self.signal=float(signal)
    self.background=float(background)
    self.constuncert=float(constuncert)/(background)
    self.lumiuncert=float(lumiuncert)/(background)
    self.lumi=float(lumi)
  def sclae_to_lumi(self,lumi):
    self.signal=(self.signal*lumi)/self.lumi
    self.background=(self.background*lumi)/self.lumi
    self.lumiuncert=self.lumiuncert*math.sqrt(self.lumi/lumi)
    self.lumi=float(lumi)
  def run_bash(self,cmd):
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    out = proc.stdout.read().strip()
    return out  #This is the stdout from the shell command
  def observation(self):
    return self.signal+self.background
  def fullunc(self):
    return math.sqrt(math.pow(self.lumiuncert,2)+math.pow(self.constuncert,2))*self.background
  def run_combine_stat(self):
    in_str="imax 1 number of channels\n"
    in_str+="jmax 1 number of background\n"
    in_str+="kmax 2 number of nuisance parameters\n"
    in_str+="observation   \t\t"+str(self.observation())+"\n"
    in_str+="bin           \t\tCn_1\t\tCn_1 \n"
    in_str+="process       \t\tSignal\t\tbkg \n"
    in_str+="process       \t\t0\t\t1 \n"
    in_str+="rate          \t\t"+str(self.signal)+"\t\t"+str(self.background)+" \n"
    in_str+="const_unc lnN \t\t-\t\t"+str(self.constuncert+1)+" \n"
    in_str+="lumi_unc lnN  \t\t-\t\t"+str(self.lumiuncert+1)+" \n"
    out_file=open("Combine.dat",'w')
    out_file.write(in_str)
    out_file.close()
    cmd_combine="combine -M ProfileLikelihood --significance Combine.dat -t -1 --expectSignal=1 2>/dev/null"
    result=self.run_bash(cmd_combine)
    #result=""
    for akt_line in result.splitlines():
      if "Significance:" in akt_line:
        return str(round(float(akt_line.split()[-1]),4))
    return "-999"
  def run_zbi_stat(self):
    n_on     = self.observation()
    mu_b_hat = self.background
    sigma_b  = self.fullunc()
    tau      = mu_b_hat / (sigma_b*sigma_b)
    n_off    = tau*mu_b_hat
    P_Bi     = TMath.BetaIncomplete(1./(1.+tau), n_on, n_off+1)
    Z_Bi     = math.sqrt(2.)*TMath.ErfInverse(1 - 2.*P_Bi)
    return str(round(Z_Bi,4))
  def run_roo_stats(self):
    n_on        = self.observation()
    mu_b_hat    = self.background
    sigma_b_rel = self.fullunc()/mu_b_hat
    signi       = RooStats.NumberCountingUtils.BinomialObsZ(n_on,mu_b_hat,sigma_b_rel)
    return str(round(signi,4))



test=SingleBin(signal,background,unc_const,unc_lumi,lumi)






print "lumi \t\tcombine\t\tzbi root \tzbi roostats\tSignal \t\tBackground \tfull_unc \tconst_unc \tlumi_unc" 
for akt_lumi in lumis:
  test.sclae_to_lumi(akt_lumi)
  print  str(akt_lumi)+"\t\t"+test.run_combine_stat()+"\t\t"+test.run_zbi_stat()+"\t\t"+test.run_roo_stats()+"\t\t"+str(round(test.signal,3))+"\t\t"+str(round(test.background,3))+"\t\t"+str(round(test.fullunc(),3))+"\t\t"+str(round(test.constuncert*test.background,3))+"\t\t"+str(round(test.lumiuncert*test.background,3))
