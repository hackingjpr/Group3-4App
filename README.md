# Idat-Shiny


- [Overview](#overview)


# Overview
This script "app.R" encodes a shiny app that upon uploading idat files, will give a relative risk score for patient samples. 
This risk score can be measured using ATRT, ECRT, or MRT (ATRT & ECRT) metagenes.

# Tutorial
## Step One
Upload your idat files (for now unzipped idat files only), including both red and green files.  
<p align="center">
<img src="https://github.com/hackingjpr/Idat-Shiny/blob/main/Tutorial/upload.png?raw=true"  
</p>
Upload a minimum of two samples. Increasing the number of samples will of course increase the length of time for the upcoming processes so we recommend ~10 sample batches. This will make looking through the results easier and will speed up the process.

## Step Two

Select your desired metagene set, this depends on wether you want your risk score to be calculated against ATRT, ECRT, or MRT metagenes.  
<p align="center">
![Metagene](https://github.com/hackingjpr/Idat-Shiny/blob/main/Tutorial/metagene.png?raw=true)  
</p>

## Step Three

Press this button:  
<p align="center">
![Generate Risk Values](https://github.com/hackingjpr/Idat-Shiny/blob/main/Tutorial/generate_risk_values.png?raw=true)  
</p>
This will start the process of generating risk values.

## Step Four: Calculating risk values

Once the calculation has been completed you should be brought to the Results tab. This tab will show a data table at the top which displays your sample names on the left and their risk vaules on the right. 
<p align="center">
![Risk Values](https://github.com/hackingjpr/Idat-Shiny/blob/main/Tutorial/Risk_values.png?raw=true)  
</p>
It will also show a graph (as shown below), the details of this graph and what it shows can be found in the corresponding paper (link)  
(insert graph example)  
There will also be a green box in the bottom right which will inform you of the currently selected sample (which will be highlighted in orange on the graph), and the selected metagene set.  
<p align="center">
![Selections](https://github.com/hackingjpr/Idat-Shiny/blob/main/Tutorial/Selections.png?raw=true)  
</p>

## Step Five: Finishing up
Once you have your results you can either reset the app by pressing the reset button:  
<p align="center">
![Reset](https://github.com/hackingjpr/Idat-Shiny/blob/main/Tutorial/Reset.png?raw=true)  
</p>
Or download your data as a CSV file (data table results), or as a PDF (data table and graph).
