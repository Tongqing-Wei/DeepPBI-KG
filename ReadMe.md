# DeepPBI-KG

### a deep learning method for prediction of phage-bacterial interactions based on key genes
#### Tips: If the ReadMe picture does not appear, make sure you are using VPN extranet access

![](https://github.com/Tongqing-Wei/DeepPBI-KG/blob/master/figure/model%20framework.png)

## Installion

##### 1. Create Conda environment

By creating a Conda environment, you can have separate environments for different projects or purposes, ensuring that the packages and dependencies within each environment do not conflict with each other.

```shell
conda create -n DeepPBI-KG python=3.9
```

##### 2. Activate your Conda environment

Once the environment is created, you can activate it and start installing packages or running programs within that environment, ensuring that the installed packages and dependencies are isolated from your system's global environment.

```shell
conda activate DeepPBI-KG
```

##### 3. Install prokka and blast using the Conda package manager

To install the prokka and blast as a dependency, once the installation is finished, you can start using prokka and blast within your Conda environment.

```shell
conda install biobuilds::prokka (conda install -c bioconda prokka)
# If the conda install not success. Download the prokka-1.12.tar.gz from https://github.com/tseemann/prokka/releases/tag/v1.12, then upload your server. Run the command:
tar zxvf prokka-1.12.tar.gz
ls prokka-1.12
cd prokka-1.12/bin

# If the prokka-1.12.tar.gz still not success, then your only download prokka 1.14 version by git.
# This causes the final predicted probability results to be slightly different from ours, as the different prokka versions make the annotated cds results slightly different.
git clone https://github.com/tseemann/prokka.git
ls prokka
cd prokka/bin

# If the installation is successful, the script below should run without any issues
prokka -h
# set database index for prokka
prokka --setupdb

conda install -c bioconda blast
blastn -h
```

##### 4. Install python modules using pip

Any of the python third-party libraries involved can be installed by pip (pip install *), where torch is installed using the following command. Refer to requirements. 

```shell
pip install torch==1.9.0+cpu torchvision==0.10.0+cpu torchaudio==0.9.0 -f https://download.pytorch.org/whl/torch_stable.html
```

## DeepPBI-KG main workflow

##### 1. Create phage database and host database for blast alignment

If you use ***git clone*** to download this project, ***all_phage_db*** and ***all_host_db*** won't get corrupted and you can skip this step (which may be slower to download).  If you download this project using ***Download Zip***, the .nsq files in ***all_phage_db*** and ***all_host_db*** may be corrupted (only a few bytes in size).  You can use the following command to create your own blast database: 

```shell
# run the python code to generate all_phage_seq.fasta and all_host_seq.fasta
python integrate_seq.py
# create database for phage
mkdir all_phage_db
makeblastdb -dbtype nucl -in all_phage_seq.fasta -parse_seqids -out ./all_phage_db/all_phage_seq.db

# create database for bacteria
mkdir all_host_db
makeblastdb -dbtype nucl -in all_host_seq.fasta -parse_seqids -out ./all_host_db/all_host_seq.db 
```

##### 2. Prokka annotation and blast alignment

Run ***prokka_blast.sh*** to get a batch of prokka annotation files. Input the fasta files of some phages (bacterium) under a folder, and output a prokka annotation folder for each phage (bacteria). Meantime, each phage (bacterium) was blast-aligned to our collected dataset to discover the phage (bacterium) with the highest homology to the query phage (bacterium). 

```shell
Usage:
  bash prokka_blast.sh [-i <FASTA_FOLDER>] [-o <OUTPUT_FOLDER>] [-p <PROKKA_PATH>] [-d <BLAST_DATABASE>] [-b <BLAST_PATH>] [-O <BLAST_OUTPUT>]
  -i FASTA_FOLDER	Input folders (some fasta files containing phages or bacteria)
  -o OUTPUT_FOLDER	Output folders (annotation result folder for each phage or bacteria generated by annotations)
  -p PROKKA_PATH	Prokka software script location path
  -d BLAST_DATABASE	Databases for blast alignment
  -b BLAST_PATH		Blastn software script location path
  -O BLAST_OUTPUT	Blastn output folders (alignment result file for each phage or bacteria)
  
Example:
  which prokka
  which blastn

  mkdir example_phage_annotation
  mkdir phage_align_result
  bash prokka_blast.sh -i ./example_phage -o ./example_phage_annotation -p ~/.conda/envs/DeepPBI-KG/bin/prokka -d ./all_phage_db/all_phage_seq.db -b ~/.conda/envs/DeepPBI-KG/bin/blastn -O ./phage_align_result
  # You should create an empty example_phage_annotation folder and phage_align_result folder in advance
  mkdir example_host_annotation
  mkdir host_align_result
  bash prokka_blast.sh -i ./example_host -o ./example_host_annotation -p ~/.conda/envs/DeepPBI-KG/bin/prokka -d ./all_host_db/all_host_seq.db -b ~/.conda/envs/DeepPBI-KG/bin/blastn -O ./host_align_result
  # You should create an empty example_host_annotation folder and host_align_result folder in advance
```

##### 3. Predictive result

Run ***DeepPBI-KG.py*** to predict PBI and directly generate prediction results (html file). Example of prediction results: 

- ***Phage*** is the phage ID to be predicted. 

- ***Bacterium*** is the bacterium ID to be predicted. 

- ***Key_gene_output*** is the prediction result of key gene model. 

- ***Wgs_output*** is the prediction result of whole genome model. 

(They are preferred to refer to the prediction result of key gene, and some key gene results are predicted to be ***nan***. Then refer to the whole genome prediction results. )

- ***Phage_key_gene_num*** is the number of key genes that a phage contains.

- ***Host_key_gene_num*** is the number of key genes that a host contains.

- ***Phage_optimum_align*** is the phage ID and E value of the optimal alignment between the query phage and the dataset collected by us, which is displayed as the hyperlink of the aligned phage ID. Clicking the hyperlink will display the specific fasta sequence of the phage. If the best alignment is empty, it is displayed as 'not hit'. 

- ***Bacterium_optimum_align*** is the bacteria ID and E value of the optimal alignment between the query bacteria and the data set collected by us, which is displayed as the hyperlink of the aligned bacteria ID. Clicking the hyperlink will display the specific fasta sequence of the bacteria. If the best alignment is empty, it is displayed as 'not hit'. 

- ***Phage_interaction_infor_in_raw_data*** is the interaction information of bacteria in the data set we collected that interacting with Phage_optimum_align. If Phage_optimum_align is empty, then interaction information is empty. 

- ***Bacterium_interaction_infor_in_raw_data*** is the interaction information of phage in the data set we collected that interacting with Bacterium_optimum_align. If Bacterium_optimum_align is empty, then interaction information is empty. 

![](https://github.com/Tongqing-Wei/DeepPBI-KG/blob/master/figure/result_example.png)

```shell
Usage:
  python DeepPBI-KG.py [--phage_annotation <PHAGE_FOLDER>] [--bacterium_annotation <BACTERIUM_FOLDER>] [--output <OUTPUT_RESULT>] [--model <MODEL_PARAMETER>] [--template <INTERMEDIATE_RESULT>] [--phage_raw_data <OUR_RAW_PHAGE_FASTA>] [--bacterium_raw_data <OUR_RAW_BACTERIUM_FASTA>] [--phage_align_res <PHAGE_ALIGN_RESULT>] [--bacterium_align_res <BACTERIUM_ALIGN_RESULT>]
  --phage_annotation PHAGE_FOLDER	Phage folder after annotation
  --bacterium_annotation BACTERIUM_FOLDER	Bacteria folder after annotation
  --output OUTPUT_RESULT	Predict result output path
  --model MODEL_PARAMETER	Related model parameter files, including model parameters and key genes list, as well as standardized scaler pkl files for input features
  --template INTERMEDIATE_RESULT	A temporary folder to store intermediate results
  --phage_raw_data OUR_RAW_PHAGE_FASTA	We collected the original dataset of phage fasta sequences
  --bacterium_raw_data OUR_RAW_BACTERIUM_FASTA	We collected the original dataset of bacterium fasta sequences
  --phage_align_res PHAGE_ALIGN_RESULT	The blast result folder for the best phage aligned with the query phage
  --bacterium_align_res BACTERIUM_ALIGN_RESULT	The blast result folder for the best bacterium aligned with the query bacterium
	
Example:
  python DeepPBI-KG.py --phage_annotation ./example_phage_annotation --bacterium_annotation ./example_host_annotation --output ./example_output --model ./model --template ./template --phage_raw_data ./phage_raw_data --bacterium_raw_data ./host_raw_data --phage_align_res ./phage_align_result --bacterium_align_res ./host_align_result
```

## Relative folder Display

- Model parameter example: 

![](https://github.com/Tongqing-Wei/DeepPBI-KG/blob/master/figure/model_parameter.png)

Where ***phage_key_gene_0.0004.csv*** and ***host_key_gene_0.0004.csv*** are the list of key genes of phage and bacteria, ***pth*** file is the model parameter, and ***pkl*** file is the standardized parameter. ***interaction_ref_infor.csv*** is all the interaction information collected for the raw data, which used to generate the Phage_interaction_infor_in_raw_data and Bacterium_interaction_infor_in_raw_data columns of  the resultant HTML files. 

- Intermediate result example: 


![](https://github.com/Tongqing-Wei/DeepPBI-KG/blob/master/figure/intermediate_result_example.png)

In ***feature_file***, feature files of phage and bacteria are displayed. ***phage_dna*** and ***phage_protein*** are fasta sequence files of DNA protein of phage cds after processing. Similarly, ***host_dna*** and ***host_protein*** are displayed.

- Output result example: 

![](https://github.com/Tongqing-Wei/DeepPBI-KG/blob/master/figure/output_result_example.png)

The ***result.csv*** and ***result.html*** is the prediction result file. The ***phage_target_align*** is the fasta sequence file of the blast alignment hit phage. The ***phage_target_interaction*** is the interaction bacteria information of the blast alignment hit phage. The ***host_target_align*** is the fasta sequence file of the blast alignment hit bacteria. The ***host_target_interaction*** is the interaction phage information of the blast alignment hit bacteria. If you want to see the complete information of the result.html, you must download the entire example_output folder. 

- Run example: 

![](https://github.com/Tongqing-Wei/DeepPBI-KG/blob/master/figure/run_example.png)

## File description



| File name                | Description                                                  |
| ------------------------ | ------------------------------------------------------------ |
| all_host_db              | sequence alignment bacteria database                         |
| all_phage_db             | sequence alignment phage database                            |
| code                     | data process and intermediate result analysis code           |
| example_host             | sample host fasta folder                                     |
| example_host_annotation  | sample host's prokka annotation result folder                |
| example_output           | example output result csv  file and html file                |
| example_phage            | sample phage fasta folder                                    |
| example_phage_annotation | sample phage's prokka annotation result folder               |
| figure                   | model illustration related figure                            |
| host_align_result        | query bacteria alignment result folder                       |
| host_raw_data            | the fasta sequence file of bacteria in our collected raw dataset |
| model                    | modle parameter related file                                 |
| phage_align_result       | query phage alignment result folder                          |
| phage_raw_data           | the fasta sequence file of phage in our collected raw dataset |
| template                 | intermediate generated feature result file                   |
| DeepPBI-KG.py            | DeepPBI-KG model predict code                                |
| integrate_seq.py         | generate all_phage_seq.fasta and all_host_seq.fasta          |
| prokka_blast.sh          | prokka annotation and sequence alignment blast code          |
| requirements.txt         | run DeepPBI-KG required related python module                |
| ReadMe.md                | ReadMe Markdown                                              |

