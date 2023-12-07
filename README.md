# HFEleStudy
For HF Electron Reconstruction Studies

Steps to run this:
1. Log in to lxplus (for the time being use lxplus7 : ```ssh username@lxplus7.cern.ch```)

2. Set up the CMSSW area:

```
$ cmsrel CMSSW_10_6_29
$ cd CMSSW_10_6_29/src
$ cmsenv
```

3. Clone the repo:
```
$ git clone git@github.com:saghosh/HFEleStudy.git
```
4. Compile:
```
$ scram b -j 8
```
5. Go tot he test area and run the code (remember to set up your voms proxy before):
```
$ cd HFEleStudy/hfe/test/
$ voms-proxy-init --voms cms
(and enter your password)
$ cmsRun ZEE_RecHit_AOD_cfg.py
```

You should get your output ROOT file to anlayse further.
