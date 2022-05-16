# Idat-Shiny


- [Overview](#overview)


# Overview
This script "app.R" encodes a shiny app that upon uploading idat files, will give a relative risk score for patient samples. 
This risk score can be measured using ATRT, ECRT, or MRT (ATRT & ECRT) metagenes.

# Tutorial
## Step One
Upload your idat files (for now unzipped idat files only), including both red and green files.  
![Upload](https://github.com/hackingjpr/Idat-Shiny/blob/main/Tutorial/upload.png?raw=true)  
Upload a minimum of two samples. Increasing the number of samples will of course increase the length of time for the upcoming processes so we recommend ~10 sample batches. This will make looking through the results easier and will speed up the process.

## Step Two

Select your desired metagene set, this depends on wether you want your risk score to be calculated against ATRT, ECRT, or MRT metagenes.  
![Metagene](https://github.com/hackingjpr/Idat-Shiny/blob/main/Tutorial/metagene.png?raw=true)

## Step Three

Press this button:  
![Generate Risk Values](https://github.com/hackingjpr/Idat-Shiny/blob/main/Tutorial/generate_risk_values.png?raw=true)  
This will start the process of generating risk values.

## Step Four

Once the calculation has been completed you should be brought to the Results tab. This tab will show a data table at the top which displays your sample names on the left and their risk vaules on the right. 
![Risk Values](https://github.com/hackingjpr/Idat-Shiny/blob/main/Tutorial/Risk_values.png?raw=true)
