# ALI-SMAD: MATLAB Space Mission Analysis and Design Tool  

ALI-SMAD is a MATLAB-based **Graphical User Interface (GUI)** for space mission analysis and design, inspired by tools like GMAT and STK.  
The software enables satellite orbit modeling, maneuver planning, and mission design with intuitive visualization and integrated analysis.  

---

## 🚀 Features  
- Orbit propagation including J2, atmospheric drag, and solar radiation pressure  
- Maneuver planning (Hohmann transfer, inclination changes, custom burns)  
- Mission goals (e.g., achieving target SMA or inclination) with ΔV calculation  
- Ground station modeling, access time analysis, and pass time prediction  
- Orbit visualization and ground track plotting  
- Integration with **STK** for result validation and comparison  

---

## 🖼️ GUI Overview  
<!-- Add figure here: Front project view with orbit visualization and ground track -->  
<img width="1366" height="722" alt="g1" src="https://github.com/user-attachments/assets/5bd75b66-aca8-4ec8-91b0-1c70dfaf9395" />

---

## 📊 Simulation Outputs  
The software provides multiple result plots and analysis options:  
- Variation of orbital elements (e.g., SMA, inclination) over time  
- Orbit trajectory visualization under perturbations (J2, drag, SRP)  

<!-- Add figure here: Example of result plots and SMA variation -->  

---<img width="1366" height="725" alt="g2" src="https://github.com/user-attachments/assets/20c4ff6e-9324-4a1d-a9d8-124dff10fcef" />


## 🔄 STK Integration  
ALI-SMAD includes the ability to export and compare results with AGI STK for validation.  

<!-- Add figure h<img width="1366" height="730" alt="g3" src="https://github.com/user-attachments/assets/fd712fde-34a8-43dc-bcf4-012c122a2c12" />
ere: Opening STK from the application and SMA variation comparison -->  

---

## 🛰️ Maneuver Planning  
The GUI supports:  
- Manual maneuver design  
- Automated Hohmann transfer calculations  
- Analysis of maneuver impact on orbital elements  

<!-- Add three figures here: Adding maneuver, example Hohmann transfer setup, and transfer results -->  
<img width="1351" height="697" alt="g4" src="https://github.com/user-attachments/assets/4c81a77c-b8e7-4c4f-b6b2-dff836a09c7a" />
<img width="1366" height="724" alt="g5" src="https://github.com/user-attachments/assets/c7f66e59-002d-475a-9872-8ebf2d91d460" />
<img width="1366" height="729" alt="g6" src="https://github.com/user-attachments/assets/91a67bff-e053-43ff-86a1-ffddacf87edf" />

---

## 🎯 Mission Goals  
Users can define target orbital elements and the tool computes required ΔV to achieve them.  
Example: Achieving an inclination change with calculated ΔV.  

<!-- Add<img width="1366" height="728" alt="g7" src="https://github.com/user-attachments/assets/5313a324-a74f-4eba-bd60-5bdb6ad35bd7" />
 figure here: Achieving mission goal with inclination change example -->  

---

## 📡 Ground Station Analysis  
Ground stations can be added with latitude, longitude, and altitude inputs.  
The software computes satellite pass times and access duration.  

<!-- Add figure here: Adding ground stations and example pass time calculation -->  
![Uploading s8.PNG…]()

---

## 📖 Reference  
- [AGI GMAT](https://gmatcentral.org/) – General Mission Analysis Tool  
- [AGI STK](https://www.agi.com/products/stk) – Systems Tool Kit  
- [VSCMG Reference – Schaub](https://hanspeterschaub.info/PapersPrivate/vscmg.pdf)  

---

## 🔑 Keywords  
Space Mission Analysis, Orbit Propagation, MATLAB GUI, STK Integration, Maneuver Planning, Ground Station Access  
