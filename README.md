The breakdown of the files is as follows: 

2txconfig.xml contains all of the required information that needs to be transmitted to the IWR1642 board to configure it for TDM-MIMO operation using both its Tx antennas.

singleframe.py is used to process signals from a single snapshot that the user triggers through TI mmWave Studio by running single_capture.lua

multitarget.py takes user inputs regarding the number of peaks that the user is trying to identify and generates a heat map along with identifying the peaks that it has found

sar_capture_only.lua is used to generate a constant stream of radar frames for the continuous data cdapture that is required for updatedsarv2.py. Note: everything from line 109 to theend is an attempt to drain the onboard RAM and reset it so the DCA1000 can keep transmitting data. However, this unfortunately did not work.

updatedsarv2.py contains the signal processing code required to take snapshots upon user trigger then combine them all into one large virtual array that is then used for backprojection SAR imaging with autofocus.
