# pcrsimulator
This downloadable app is able to simulate polymerase chain reactions in DNA. It was created in order to ease synthetic biologist's lives when designing primers and conducting an experiment in silico. I created this platform for Dr. Chris Anderson's Synthetic Biology Lab course at UC Berkeley - BioE 140L.

Link to download the app: https://drive.google.com/file/d/1IuvCjOeDEQXbEGdjWCy6XAOluWbzIipS/view?usp=sharing

**Functionalities**

This app takes three inputs: a forward primer, a reverse primer, and a template strand. It ensures that the three oligonucleotides only contain DNA or IUPAC codes, and checks to see if the primers have a high enough melting temperature in order to not degrade during the PCR process. Then, it finds the annealing region on the 3' end of the primer by starting with the first six bases and continuing in the 5' direction until it does not match the template strand anymore. The number six was chosen here because that is theoretically the minimum number of bases you need to anneal, though a much larger annealing region is more favorable. The PCR Simulator then includes additional bases on the 5' end of the primers if necessary. It treats the template DNA as double stranded and circular. If someone would like to do iPCR (inverse PCR), they can. Looking at the image below, one would simply insert the 'REV' primer as the Forward Primer, and the 'FW' primer as the reverse primer. The template strand could be inserted as is.

![Figure depicting process of iPCR.](http://atlasgeneticsoncology.org/Deep/Images/LDI-PCRinCancerFig1.png)
http://atlasgeneticsoncology.org/Deep/LDI-PCRinCancerID20087.html

**Notes**

The original assignment was to create an Excel Add-in that functioned as a PCR simulator. When looking into this, I found that UDFs (user-defined functions) were prohibited in Excel 2016 for Mac. Thus, the only way to create functions on a Mac would be to code them in visual basic in the VBE (visual basic editor) that is built into Excel. If a future student attempts to solve this, they would need a Windows or Linux machine, and I would advise them to use the downloadable python module xlwings, which allows them to code in python and simply point to their python file within Excel's VBE. Since I was unable to create an Excel add-in, I simply created a Macintosh app that has the same functionality. If someone with Windows wanted to use it, I uploaded my python code to this GitHub repository and they could simply download that and wrap it in Excel using xlwings. 

**Acknowledgements**

I used PyCharm to write and debug my source code, and PyQt5 to create a GUI. I used fbs to build the installer for the app. I would like to thank Professor Anderson for his guidance along the way.
