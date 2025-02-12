# ShinyIBL

This is the landing page for ShinyIBL, a Shiny app that allows a user to run basic decisions-from-experience simulations based on the Instance-Based Learning Theory (Gonzalez et al, 2003). 

## Demo 

The app can be demoed at http://ddmlab.com:3838/shinyibl/.

## Run Remotely

### Install R and RStudio

To install R and RStudio, follow the instructions here for your operating system: https://rstudio-education.github.io/hopr/starting.html

It is important to install R first, then RStudio.

### Run App

1. After installing both R and RStudio, open the RStudio program.
2. Download the code for this repository by clicking "Code" and downloading as a ZIP. Alternatively, you can clone the repository if you are comfortable with GitHub.
3. Unzip the folder. You will now work with the unzipped folder.
4. Open RStudio. 
5. To open the app and run it, there are two options. 

- Option 1 (Recommended): Create a new project by clicking "File" -> "New Project" -> "Existing Directory" and selecting the unzipped folder. This will set the working directory to the unzipped folder.
- Option 2: Set the working directory to be the unzipped folder. There are two main ways to do this. 
    -  Run the following in the console, substituting YOURFILEPATHHERE with the file path to the folder. Make sure to include the quotation marks, and also note that R uses forward slashes (/) as opposed to back slashes (\).
```
setwd("YOURFILEPATHHERE")
```
If you saved the folder to your downloads folder, for example, you would use something along the lines of the following:
```
setwd("/Users/YOURUSERNAME/Downloads/shinyIBL-main")
```
    - Or, navigate to the folder using the bottom right panel of RStudio by clicking through the "Files" tab. Then, click "More" and then "Set As Working Directory". You can check the current working directory by running the following in the console:
```
getwd()
```
6. Once you have created a project or set the working directory, open the app.R file in RStudio.
7. Run the following in the console and follow the instructions to install the necessary packages:
```
install.packages("shiny")
install.packages("ggvis")
install.packages("dplyr")
install.packages("reshape2")
install.packages("tidyr")
install.packages("rhandsontable")
install.packages("data.table")
```
8. Type `shiny::runApp()` in the console or click the "Run" button on the upper right. The app should run locally.

Email Erin Bugbee (ebugbee@andrew.cmu.edu) with any questions about the app.