import streamlit as st #allows for building web application
import matplotlib.pyplot as plt #allows for plotting
import pandas as pd #allows for data handling
import numpy as np #allows for numerical operations
import spectrum_utils.plot as sup
import spectrum_utils.spectrum as sus
import plotly.tools as tls

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
    plt.figure(figsize=(10, 6))
    plt.stem(mzs, intensities, basefmt=" ", use_line_collection=True)
    plt.title(f"MS2 Spectrum for Scan {scanNum}")
    plt.xlabel("m/z", fontsize=11)
    plt.ylabel("Intensity", fontsize=11)
    st.pyplot(plt)

@st.cache_data
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

        df = pd.DataFrame(scans, columns = ["Scan Number", "Spectrum ID", "PEPMASS Number", "Charge State", "SMILES ID"])
        
        # Dropdown menu for selecting scan number
        scan_input = st.selectbox("Select Scan Number to view MS2 Spectrum", options=scan_nums)
        
        # Display the DataFrame in an expander to make it more compressed
        with st.expander("Show Scans"):
            st.dataframe(df) #Outputs Dataframe of Scans

        if scan_input:
            user_scan = next(scan for scan in scans if scan["Scan Number"] == scan_input)
            mz_filtered, sqrt_filtered, normal_filtered = peak_filtering(scan_input)
            if st.button("View Unfiltered Spectrum"):
                peak_visual(user_scan["m/z data"], user_scan["intensity data"], str(user_scan["Scan Number"]), user_scan["PEPMASS Number"], user_scan["Charge State"])
            if st.button("View Filtered Spectrum - Normalized"):
                peak_visual(mz_filtered, normal_filtered, str(user_scan["Scan Number"]), user_scan["PEPMASS Number"], user_scan["Charge State"])
            if st.button("View Filtered Spectrum - Square Root Normalized"):
                peak_visual(mz_filtered, sqrt_filtered, str(user_scan["Scan Number"]), user_scan["PEPMASS Number"], user_scan["Charge State"])
