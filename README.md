## TauHLT

```
cmsrel CMSSW_12_1_0
cd CMSSW_12_1_0/src  
cmsenv    
git cms-addpkg HLTrigger/Configuration    
git clone git@github.com:sarafiorendi/TauHLT.git
cd TauHLT/  
git checkout main   
cd ..  
scramv1 b   
```  

To produce ntuples:
```   
cmsRun hltTauNtuples.py
```

This gist customize the hlt configuration to re-run the hlt and the ntuplizer displaced tau studies
(update link to gist)
https://gist.github.com/sarafiorendi/02b4a43c28766ec7beddabe3f5ae036e
  