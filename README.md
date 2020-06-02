# HGCNose

Author: K.Yoon

## Setup
### 1. CMSSW working directory
If already done, skip to 2.

### 2. Environment variables
In bash, open config.sh and edit `DIR_DATA` to whichever directory your simulated data is stored.
```bash
. config.sh
```

## Analysis
Here is the list of analysis modules.
* EnergyResolution | Returns energy resolution for a given PID and eta range. Outputs plots in output/test.root (because it is being tested).
