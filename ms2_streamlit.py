import streamlit as st #allows for building web application
import matplotlib.pyplot as plt #allows for plotting
import spectrum_utils.plot as sup
import spectrum_utils.spectrum as sus
import plotly.tools as tls
import pandas as pd #allows for data handling
import numpy as np #allows for numerical operations

def peak_filtering(user_scan): #DONE
    mz_array = user_scan["m/z data"]
    intensity_array = user_scan["intensity data"]
    filtered_mz = []
    filtered_intensities = []
    pepmass_val = float(user_scan["PEPMASS Number"])

        #basic peak filtering
    for i, mz in enumerate(mz_array):
        peak_range = [j for j in range(len(mz_array)) if abs(mz_array[j] - mz) <= 25]
        sorted_range = sorted(peak_range, key=lambda j: intensity_array[j], reverse=True)
        if i in sorted_range[:6]:
            if abs(mz - pepmass_val) > 17:
                filtered_mz.append(mz)
                filtered_intensities.append(intensity_array[i])

    sqrt_data, normalized_data = peak_normalizing(filtered_intensities)
        
    return filtered_mz, sqrt_data, normalized_data
    
def peak_normalizing(filtered_intensities): #DONE

    normalized_intensities = np.copy(filtered_intensities) / np.linalg.norm(filtered_intensities)
    sqrt_intensities = np.sqrt(normalized_intensities)

    return sqrt_intensities, normalized_intensities

def peak_visual(mzs, intensities, scanNum, pepmass, charge):

    spectrum = sus.MsmsSpectrum(mz = mzs, intensity=intensities, identifier=scanNum, precursor_mz=pepmass, precursor_charge=charge)
    sup.spectrum(spectrum)
    plt.title(f"MS2 Spectrum for Scan {scanNum}")
    plt.xlabel("m/z", fontsize = 11)
    plt.ylabel("Intensity", fontsize = 11)
    fig = plt.gcf()
    plotly_fig = tls.mpl_to_plotly(fig)
    plotly_fig.update_traces(hoverinfo="x+y")
    
    st.plotly_chart(plotly_fig)  

def read_mgf_file(mgf_file): #DONE
    scans = [] #list to store parsed scans
    current_scan = None
    scan_numbers = []

    file = mgf_file.read().decode('utf-8').splitlines()
    for line in file: #for each line in the file
        line = line.strip()
        if line == "BEGIN IONS": #beginning of scan 
            current_scan = {"Scan Number": 0, "Spectrum ID": '', "PEPMASS Number": 0.0, "Charge State": 0, "SMILES ID": '', "peaks": [], "m/z data": [], "intensity data": []} #initializes new scan with keys
        elif line == "END IONS": #end of scan
            if current_scan:
                scans.append(current_scan) #adding current scan to total scans
                current_scan = None #Reseting for next scan in mgf file
        elif current_scan is not None: #if scan has begun
            if "=" in line:
                data = line.split('=', 1) #limits line split to first '='
                if len(data) == 2:
                    key, value = data
                    if key == "SCANS":
                        current_scan["Scan Number"] = int(value)
                        scan_numbers.append(int(value))
                    elif key == "SPECTRUMID":
                        current_scan['Spectrum ID'] = str(value)
                    elif key == "PEPMASS":
                        current_scan["PEPMASS Number"] = float(value)
                    elif key == "CHARGE":
                        current_scan["Charge State"] = int(value)
                    elif key == 'SMILES':
                        current_scan["SMILES ID"] = str(value.strip())
                    else: 
                        continue
            else: #must be peak data
                try: 
                    data2 = line.split()
                    if len(data2) == 2:
                        mz, intensity = data2
                        current_scan["peaks"].append((float(mz), float(intensity)))
                        current_scan["m/z data"].append(float(mz))
                        current_scan["intensity data"].append(float(intensity))
                except ValueError: 
                    print(f"Skipping unreadable data in line: '{line}")
                    continue
        
    return scans, scan_numbers
    
if __name__ == "__main__":
    st.title("MS2 Scan") #App title

    mgf_file = st.file_uploader("Choose a file", type="mgf") #allows user to upload file
    
    if mgf_file is not None: #ensures that the file holds a valid reference
        scans, scan_nums = read_mgf_file(mgf_file)
        #Search Box For Scan Number -> Searches for Searched Scan Number within file
        search_scan_input = st.text_input("Enter Scan Number to view MS2 Spectrum")
        user_scan = {}
        if search_scan_input:
            if search_scan_input.isdigit():
                search_scan = int(search_scan_input)
                if search_scan in scan_nums:
                    user_scan = next(scan for scan in scans if scan["Scan Number"] == search_scan)
                    st.write("Scan Number Found Within File - Use Filters to View Data")
                else:
                    st.write("Scan Number Not Found. Please enter a Valid Scan Number")
                    search_scan = None
            else:
                st.write("Scan Number Not Found. Please enter a Valid Scan Number")
                search_scan = None
        else: 
            st.write("(Please enter Scan Number)")
        # Output all Scan Numbers using a Checkbox
        if st.checkbox("Show Scan Numbers"):
            # st.write("Scan Numbers:")
            if scan_nums is not None:
                st.write(f"{scan_nums}")
        if user_scan:
            #Outputs MS2 Meta Data
            meta_data = False
            if st.checkbox("View Meta-Data"):
                meta_data = True
            #Allows for Peak Filtering     
            mz_filtered, sqrt_filtered, normal_filtered = peak_filtering(user_scan)

            #Normalizes Peaks using a Checkbox
            filtered = False
            sqrt_data = False

            if st.checkbox("Filter Peaks (Normal)"):
                filtered = True
            if st.checkbox("Filter Peaks (Square Root)"):
                filtered = True
                sqrt_data = True  

            #Visualizes MS2 Spectrum
            if st.button("Visualize MS2 Spectrum"):
                if meta_data:
                    st.write("Scan Number: ", user_scan["Scan Number"])
                    st.write("Spectrum ID: ", user_scan["Spectrum ID"])
                    st.write("PEPMASS Number: ", user_scan["PEPMASS Number"])
                    st.write("Charge State: ", user_scan["Charge State"])
                    st.write("SMILES ID: ", user_scan["SMILES ID"])
                else: 
                    st.write("Meta Data Not Shown")   
                if filtered:
                    if sqrt_data:
                        peak_visual(mz_filtered, sqrt_filtered, str(user_scan["Scan Number"]), user_scan["PEPMASS Number"], user_scan["Charge State"])
                    else:
                        peak_visual(mz_filtered, normal_filtered, str(user_scan["Scan Number"]), user_scan["PEPMASS Number"], user_scan["Charge State"])
                else: #gives unfiltered plot    
                    peak_visual(user_scan["m/z data"], user_scan["intensity data"], str(user_scan["Scan Number"]), user_scan["PEPMASS Number"], user_scan["Charge State"])     
                
                
                    





    
        
         

    
