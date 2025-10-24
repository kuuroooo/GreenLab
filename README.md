

### Description of Files and Structure 
- The jobRunner.sh & RunnerConfig.py, which were used to run the experiment, are to be found in the benchmarks directory

- Run table can be found in experiments/benchmarks_energy_analysis/

- The script used for analysis can be found under analysis.R

***
### Target Machine Setup


Enable SSH on the target system, (MacBook)
#### Through system settings:
1. Open **System Settings → General → Sharing**
2. Find and toggle **Remote Login**
3. When it turns green (enabled), it automatically starts the SSH server

#### Copy headless.py
From the Experiment Orchestrator, copy the *headless.py* file into the root directory of your target device. 
This can be done using scp if you have already successfully setup both devices SSH connection.

```
scp path/to/headless.py target_user@target_ip/
```
For simplicity we will leave the *headless.py* in the root directory

Now place all repository files onto the target system, install energi bridge and all dependencies (see experiment-runner), once this is done your setup for the target system is complete.

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
Now that both SSH and the target system are setup, you will need to adjust the relevant file paths in a variety of scripts. 


***
## Running the experiment



***
## Analyzing the data
Given that you have installed R, if not do so now, on the system where you have your EXPERIMENTS and Results folder you should run the analysis.R file in the following way:
```
Rscript analysis.R
```
This should install all necessary R dependencies and generate the graphs, as long as it is located in the same directory as EXPERIMENTS and Results

Congratulations you are done !


