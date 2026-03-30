# Solving Trade Models for Policy Simulations

This project provides tools and guides for solving trade models for researchers, economists and students.

Human generated solution code is available in Python, Julia, Matlab, and R for Eaton&Kortum(2002)-like trade models. Run the script "AO_Extended.x" for the full model and "EK_Simple.x" for the basic model. Each script is self-contained and just requires an included data file. The solution algorithm is linear and easy to follow. Processed data, quick reference for equations, sector list and country list are included in the folder. See the instructions below.

---

## Run a policy scenario online: Trade Policy Simulator 

If you are in a rush or don't want to deal with programming, dashboards are available to solve the full model in real time on the server. Plug a tariff change or a trade cost scenario to see the impact on real disposable income, GDP, exports and imports globally. Different from many other dashboards, the model is solved on-demand so any scenario is possible (subject to interface constraints).

***Trade Policy Simulator: Tariffs***\
https://www.worldbank.org/en/data/interactive/2025/12/15/trade-simulator-tariffs

***Trade Policy Simulator: Trade Costs***\
https://www.worldbank.org/en/data/interactive/2025/12/15/trade-simulator-trade-costs

---

***Paper***:\
Artuc and Ortega (2026). "International Trade Policy and Quantitative Models: A Practitioner’s Guide," World Bank PRWP.\
(forthcoming, the link will be available soon)

***Full replication package for the paper***:\
https://reproducibility.worldbank.org/catalog/511

Please note that we might update the code from time to time. The code here and the generated results might be slightly different than the replication package and the paper, due to small changes in data processing assumptions, approximation differences, or improvements for readibility/efficiency.

***Data sources:***\
Base year 2022 with 21 sectors and 81 countries.\
Trade data: TiVA-ICIO 2025, January 2026 update.\
Tariff data: Constantinescu (2026), November 2025 update.


***Links***\
Erhan: http://www.artuc.org  
Johan: https://www.linkedin.com/in/johan-ortega-h-829665123/

---

### Instruction for offline solution

Run one of the scripts below. Equations from the model are clearly marked within the code including 
references to their numbers in the paper and the corresponding steps in the algorithm. 
Extended model requires “data_tiva25_tariff.mat” file, 
and the simple model requires “data_tiva25_simple.mat” file. Both data files are provided. Results will be saved as a csv file showing changes in disposable income (i.e. income excluding tariff revenue), total output, exports and imports, all in real terms.


| Script | Language |Model | 
| :--- | :--- | :--- |
| AO_Extended.py |Python |Extended |
| AO_Extended.jl    |Julia |Extended |
| AO_Extended.m | Matlab |Extended |
| AO_Extended.R    |R    |Extended |
| EK_Simple.py |Python |EK |
| EK_Simple.jl    |Julia |EK |
| EK_Simple.m |Matlab |EK |
| EK_Simple.R    |R    |EK |



