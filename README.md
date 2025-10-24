

### Description of Files and Structure 
- The jobRunner.sh & RunnerConfig.py, which were used to run the experiment, are to be found in the targetSystem/benchmarks directory

- Run table can be found in targetSystem/experiments/benchmarks_energy_analysis/

- The script used for analysis can be found under targetSystem/analysis.R

- The jobRunner.sh for the ExperimentOrchestrator can be found under ExperimentOrchestrator/jobRunner.sh

- The WarmupScript.py used to pre heat the target system, RunnerConfig.py used to remotely execute the experiment runner can all be found under the ExperimentOrchestrator/

- The headless.sh used to reduce the energy consumption of the target system by shutting down unnecessary applications can be found under: targetSystem/headless.sh
***
### Target Machine Setup


## 1.)Enable SSH on the target system, (MacBook)
#### Through system settings:
1. Open **System Settings → General → Sharing**
2. Find and toggle **Remote Login**
3. When it turns green (enabled), it automatically starts the SSH server

#### 2.)Copy headless.sh
From the targetSystem, copy the *headless.sh* file into the root directory of your target device. 

```
mv targetSystem/headless.sh /
```
For simplicity we will leave the *headless.sh* in the root directory


Most likely you will need to make the *headless.sh* file executable using:
```
sudo chmod +x headless.sh
```
On your target device

Now place all repository files onto the target system, install energi bridge and all dependencies (see experiment-runner), once this is done your setup for the target system is complete.

### Install the Dependencies
After successfully setting up the experiment runner, *see experiment runner*
You should have installed the majority of requirements needed.

***
### Experiment Orchestrator Setup

#### Steps

### 1.) Generate an SSH Key Pair on your Experiment Orchestrator with your target machine (server) where the experiments will run on

Open a terminal on your experiment orchestrator and run: 

```bash
ssh-keygen -t ed25519 -C "target_username@target_ip"
````

- Save location: press Enter to accept default (~/.ssh/id_ed25519)
- Passphrase: press Enter for no passphrase (or add one for extra security)

This creates two files:
- ~/.ssh/id_ed25519        # Private key (keep secret)
- ~/.ssh/id_ed25519.pub    # Public key (shareable)

Copy the key to your target machine (In our case the MacBook) 
```
ssh-copy-id user@ip
```
Now test the connection to the target machine using:
```
ssh user@ip
```
If everything works we now have to create the ssh shortcut beginning with:
Edit (or create) this file on your experiment orchestrator
```
nano ~/.ssh/config
```
Add:
```
Host macbook
  Hostname TARGET IP GOES HERE
  User YOUR USER GOES HERE
  IdentityFile ~/.ssh/id_ed25519
```
Now we can test the shortcut:
```
ssh macbook
```
If successful, the SSH configuration for the experiment orchestrator should now be setup properly and scripts can now use this shortcut.

### 2.) Adjust path variables
Now that both SSH and the target system are setup, you will need to adjust the relevant file paths in the .env file for the scripts to work.

Navigate to your Experiment Orchestrators directory and edit or create a .env file, you can do this through:
```
touch .env (if .env doenst exist)
nano .env
```
In the *.env* file add the following two lines:
```
SSH_COMMAND=ssh macbook
REMOTE_REPO=/Users/kellywang/Documents/compSci/p1/greenLab/GreenLab
```
Depending on how you followed the ssh setup in step 1, make sure to use the appropriate shortcut name you setup. 

Furthermore make sure to replace the path for the REMOTE_REPO with the relevant REMOTE_REPO path, where your experiment is located on 
the target machine respectively

Now you are good to go !

### 3.) Install any required dependencies
On the Experiment orchestrator make sure to run:
```
pip install -r requirements.txt
```
This should cover most of the requirements needed...

***
## Running the experiment
Once you have followed the previous steps navigate to the Experiment Orchestrator directory on your experiment orchestrator device.
(Optional) Within the jobRunner.sh script you can set the number of experiment runs, default is 10

In order to run the experiment, go to the Experiment Orchestrator Folder and simply execute the jobRunner.sh script using:
```
./jobRunner.sh
```
If you receive a permission error you might need to make the script executable using 
```
sudo chmod +x jobRunner.sh
```

Done, The target PC should now warm up and begin running the experiments, this can take up to several hours, depending on your experiment number setting, 
default is 10 (Used within this report)


Once finished you can now retrieve the data from the Experiment Orchestrator, it will be found in a variety of folders.

***
## Analyzing the data
Given that you have installed R, if not do so now, on the system where you have your EXPERIMENTS and Results folder you should run the analysis.R file in the following way:
```
Rscript analysis.R
```
This should install all necessary R dependencies and generate the graphs, as long as it is located in the same directory as EXPERIMENTS and Results

Congratulations you are done !


