
import matplotlib.pyplot as plt
import IslandSystem as IS
import pandas as pd

# ------------------------- Läser all data --------------------------
# File paths for the two Excel files
file_paths = ['elkonsumtion.xlsx', 'power_output_solceller_hustak.xlsx', 'power_output_solceller_mark.xlsx', 'energy_production_wind2.xlsx', 'average_wind_speeds.xlsx']

# Specify the columns you want to read from each file
columns_to_read = {
    'elkonsumtion': ['Unnamed: 1', 'Avläst elanvändning per timme'],
    'power_output_solceller_hustak': ['Date and time', 'Cell temperature', 'Module efficiency', 'PV system output'],
    'power_output_solceller_mark': ['Date and time', 'Cell temperature', 'Module efficiency', 'PV system output'],
    'energy_production_wind2': ['Date', 'Power Output (W)'],
    'average_wind_speeds': ['Month-Day-Time', 'Average Wind Speed (m/s)']
}

# Create an empty dictionary to store dataframes
dfs = {}

# Read each Excel file into a dataframe and store it in the dictionary
for file_path in file_paths:
    # Extract filename without extension
    filename = file_path.split('.')[0]

    # Read Excel file into a dataframe using the openpyxl engine
    df = pd.read_excel(file_path, engine='openpyxl', usecols=columns_to_read.get(filename))

    # Store the dataframe in the dictionary with filename as the key
    dfs[filename] = df

#lägg till Aeosol
file_path_A = 'energy_production_wind2.xlsx'
filename_A = file_path_A.split('.')[0]+'_A'
dfs[filename_A] = pd.read_excel(file_path_A, engine='openpyxl', usecols=['Date', 'Power Output (Wh)'], sheet_name='Aeolos_V-3000')

# -------- gör om från halvtimme till timme vind -----------

x = 0
energy_production_wind_hour = []
while x < 17520: # byt denna till 17520
    column_index = dfs['energy_production_wind2'].columns.get_loc('Power Output (W)')
    datapoint1 = dfs['energy_production_wind2'].iloc[x, column_index]
    datapoint2 = dfs['energy_production_wind2'].iloc[x+1, column_index]
    energy_production_wind_hour.append(datapoint1 + datapoint2)
    x += 2

x = 0
energy_production_wind_hour_A = []
while x < 17520: # byt denna till 17520
    column_index = dfs['energy_production_wind2_A'].columns.get_loc('Power Output (Wh)')
    datapoint1 = dfs['energy_production_wind2_A'].iloc[x, column_index]
    datapoint2 = dfs['energy_production_wind2_A'].iloc[x+1, column_index]
    energy_production_wind_hour_A.append(datapoint1 + datapoint2)
    x += 2

# Initialize the sum_power list with zeros
sum_power_tmvW = [] # tak, mark, vind (Windstar)
sum_power_tmvA = [] # tak, mark, vind (Aeolos)
sum_power_tm = [] # tak, mark
sum_power_mvW = [] # mark, vind (Aeolos)
sum_power_mvA = []
sum_power_tvW = [] # tak, vind (W)
sum_power_tvA = []
sum_power_t = [] # tak
sum_power_m = [] # mark
sum_power_vW = [] # vind Windstar
sum_power_vA = [] # vind Aeolos


# Loop through each index and calculate sum_power
for i in range(8760):
    sum_power_tmvW.append(dfs['power_output_solceller_hustak'].iloc[i, dfs['power_output_solceller_hustak'].columns.get_loc('PV system output')]/1000 + dfs['power_output_solceller_mark'].iloc[i, dfs['power_output_solceller_mark'].columns.get_loc('PV system output')]/1000 + energy_production_wind_hour[i]/1000 - dfs['elkonsumtion'].iloc[i, dfs['elkonsumtion'].columns.get_loc('Unnamed: 1')])

    sum_power_tmvA.append(dfs['power_output_solceller_hustak'].iloc[i, dfs['power_output_solceller_hustak'].columns.get_loc('PV system output')]/1000 + dfs['power_output_solceller_mark'].iloc[i, dfs['power_output_solceller_mark'].columns.get_loc('PV system output')]/1000 + energy_production_wind_hour_A[i]/1000 - dfs['elkonsumtion'].iloc[i, dfs['elkonsumtion'].columns.get_loc('Unnamed: 1')])

    sum_power_tm.append(dfs['power_output_solceller_hustak'].iloc[i, dfs['power_output_solceller_hustak'].columns.get_loc('PV system output')]/1000 + dfs['power_output_solceller_mark'].iloc[i, dfs['power_output_solceller_mark'].columns.get_loc('PV system output')]/1000 - dfs['elkonsumtion'].iloc[i, dfs['elkonsumtion'].columns.get_loc('Unnamed: 1')])

    # För solenergi bara från mark med vindkraft
    sum_power_mvW.append(dfs['power_output_solceller_mark'].iloc[i, dfs['power_output_solceller_mark'].columns.get_loc('PV system output')]/1000 + energy_production_wind_hour[i]/1000 - dfs['elkonsumtion'].iloc[i, dfs['elkonsumtion'].columns.get_loc('Unnamed: 1')])

    sum_power_mvA.append(dfs['power_output_solceller_mark'].iloc[i, dfs['power_output_solceller_mark'].columns.get_loc('PV system output')]/1000 + energy_production_wind_hour_A[i]/1000 - dfs['elkonsumtion'].iloc[i, dfs['elkonsumtion'].columns.get_loc('Unnamed: 1')])

    sum_power_tvW.append(dfs['power_output_solceller_hustak'].iloc[i, dfs['power_output_solceller_hustak'].columns.get_loc('PV system output')]/1000 + energy_production_wind_hour[i]/1000 - dfs['elkonsumtion'].iloc[i, dfs['elkonsumtion'].columns.get_loc('Unnamed: 1')])

    sum_power_tvA.append(dfs['power_output_solceller_hustak'].iloc[i, dfs['power_output_solceller_hustak'].columns.get_loc('PV system output')]/1000 + energy_production_wind_hour_A[i]/1000 - dfs['elkonsumtion'].iloc[i, dfs['elkonsumtion'].columns.get_loc('Unnamed: 1')])

    sum_power_t.append(dfs['power_output_solceller_hustak'].iloc[i, dfs['power_output_solceller_hustak'].columns.get_loc('PV system output')]/1000 - dfs['elkonsumtion'].iloc[i, dfs['elkonsumtion'].columns.get_loc('Unnamed: 1')])

    sum_power_m.append(dfs['power_output_solceller_mark'].iloc[i, dfs['power_output_solceller_mark'].columns.get_loc('PV system output')]/1000 - dfs['elkonsumtion'].iloc[i, dfs['elkonsumtion'].columns.get_loc('Unnamed: 1')])

    sum_power_vW.append(energy_production_wind_hour[i] / 1000 - dfs['elkonsumtion'].iloc[i, dfs['elkonsumtion'].columns.get_loc('Unnamed: 1')])

    sum_power_vA.append(energy_production_wind_hour_A[i] / 1000 - dfs['elkonsumtion'].iloc[i, dfs['elkonsumtion'].columns.get_loc('Unnamed: 1')])

sum_power_list = [sum_power_tmvW, sum_power_tmvA, sum_power_tm, sum_power_mvW, sum_power_mvA, sum_power_tvW, sum_power_tvA, sum_power_t, sum_power_m, sum_power_vW, sum_power_vA]
sum_power_list_names = ['sum_power_tmvW', 'sum_power_tmvA', 'sum_power_tm', 'sum_power_mvW', 'sum_power_mvA', 'sum_power_tvW', 'sum_power_tvA', 'sum_power_t', 'sum_power_m', 'sum_power_vW', 'sum_power_vA']

results = {}
batteri_quant = [0, 1, 2, 3, 4, 5]
H_lagring_quant = [0, 1]

# ----------- Loopar så att det blir fler år
years = 5  # Antal år simuleras. Bra för att se laddning av vätgas
power_step = 0
for power in sum_power_list:
    sum_power = IS.data_extender(power, years)
    # ---------------------------------------------
    # vind_quant = []
    name = sum_power_list_names[power_step]
    power_step += 1
    for x in batteri_quant:
        for y in H_lagring_quant:
            batteri, H_lagring, buy_electricity, sell_electricity = IS.system(sum_power, x, y)
            result_key = f"result_{name}_{x}_{y}"  # Dynamically generate result key
            results[result_key] = (batteri, H_lagring, buy_electricity, sell_electricity)


summary = {}
# Iterate over the keys in the results dictionary
for key, value in results.items():
    total_buy_electricity_sum = 0
    total_sell_electricity_sum = 0

    # Extract the 'batteri' values from the current result
    batteri_values = value[0] if isinstance(value[0], list) else [value[0]]
    H_lagring_values = value[1] if isinstance(value[1], list) else [value[1]]
    buy_electricity_values = value[2] if isinstance(value[2], list) else [value[2]]
    sell_electricity_values = value[3] if isinstance(value[3], list) else [value[3]]

    mean_batteri = sum(batteri_values) / len(batteri_values)
    mean_H_lagring = sum(H_lagring_values) / len(H_lagring_values)
    total_buy_electricity_sum += sum(buy_electricity_values)
    total_sell_electricity_sum += sum(sell_electricity_values)

    value_list = [mean_batteri, mean_H_lagring, total_buy_electricity_sum, total_sell_electricity_sum]
    summary_key = f"summary_{key}"  # Dynamically generate result key
    summary[summary_key] = value_list

#---------------------------------------------------------------------------------------------

tmvW_matrix_batteri = [['tmvW',0,1,2,3,4,5], [0,0,0,0,0,0,0], [1,0,0,0,0,0,0]]
tmvW_matrix_H_lagring = [['tmvW',0,1,2,3,4,5], [0,0,0,0,0,0,0], [1,0,0,0,0,0,0]]
tmvW_matrix_buy = [['tmvW',0,1,2,3,4,5], [0,0,0,0,0,0,0], [1,0,0,0,0,0,0]]
tmvW_matrix_sell = [['tmvW',0,1,2,3,4,5], [0,0,0,0,0,0,0], [1,0,0,0,0,0,0]]

tmvA_matrix_batteri = [['tmvA',0,1,2,3,4,5], [0,0,0,0,0,0,0], [1,0,0,0,0,0,0]]
tmvA_matrix_H_lagring = [['tmvA',0,1,2,3,4,5], [0,0,0,0,0,0,0], [1,0,0,0,0,0,0]]
tmvA_matrix_buy = [['tmvA',0,1,2,3,4,5], [0,0,0,0,0,0,0], [1,0,0,0,0,0,0]]
tmvA_matrix_sell = [['tmvA',0,1,2,3,4,5], [0,0,0,0,0,0,0], [1,0,0,0,0,0,0]]

tm_matrix_batteri = [['tm',0,1,2,3,4,5], [0,0,0,0,0,0,0], [1,0,0,0,0,0,0]]
tm_matrix_H_lagring = [['tm',0,1,2,3,4,5], [0,0,0,0,0,0,0], [1,0,0,0,0,0,0]]
tm_matrix_buy = [['tm',0,1,2,3,4,5], [0,0,0,0,0,0,0], [1,0,0,0,0,0,0]]
tm_matrix_sell = [['tm',0,1,2,3,4,5], [0,0,0,0,0,0,0], [1,0,0,0,0,0,0]]

mvW_matrix_batteri = [['mvW',0,1,2,3,4,5], [0,0,0,0,0,0,0], [1,0,0,0,0,0,0]]
mvW_matrix_H_lagring = [['mvW',0,1,2,3,4,5], [0,0,0,0,0,0,0], [1,0,0,0,0,0,0]]
mvW_matrix_buy = [['mvW',0,1,2,3,4,5], [0,0,0,0,0,0,0], [1,0,0,0,0,0,0]]
mvW_matrix_sell = [['mvW',0,1,2,3,4,5], [0,0,0,0,0,0,0], [1,0,0,0,0,0,0]]

mvA_matrix_batteri = [['mvA',0,1,2,3,4,5], [0,0,0,0,0,0,0], [1,0,0,0,0,0,0]]
mvA_matrix_H_lagring = [['mvA',0,1,2,3,4,5], [0,0,0,0,0,0,0], [1,0,0,0,0,0,0]]
mvA_matrix_buy = [['mvA',0,1,2,3,4,5], [0,0,0,0,0,0,0], [1,0,0,0,0,0,0]]
mvA_matrix_sell = [['mvA',0,1,2,3,4,5], [0,0,0,0,0,0,0], [1,0,0,0,0,0,0]]

tvW_matrix_batteri = [['tvW',0,1,2,3,4,5], [0,0,0,0,0,0,0], [1,0,0,0,0,0,0]]
tvW_matrix_H_lagring = [['tvW',0,1,2,3,4,5], [0,0,0,0,0,0,0], [1,0,0,0,0,0,0]]
tvW_matrix_buy = [['tvW',0,1,2,3,4,5], [0,0,0,0,0,0,0], [1,0,0,0,0,0,0]]
tvW_matrix_sell = [['tvW',0,1,2,3,4,5], [0,0,0,0,0,0,0], [1,0,0,0,0,0,0]]

tvA_matrix_batteri = [['tvA',0,1,2,3,4,5], [0,0,0,0,0,0,0], [1,0,0,0,0,0,0]]
tvA_matrix_H_lagring = [['tvA',0,1,2,3,4,5], [0,0,0,0,0,0,0], [1,0,0,0,0,0,0]]
tvA_matrix_buy = [['tvA',0,1,2,3,4,5], [0,0,0,0,0,0,0], [1,0,0,0,0,0,0]]
tvA_matrix_sell = [['tvA',0,1,2,3,4,5], [0,0,0,0,0,0,0], [1,0,0,0,0,0,0]]

t_matrix_batteri = [['t',0,1,2,3,4,5], [0,0,0,0,0,0,0], [1,0,0,0,0,0,0]]
t_matrix_H_lagring = [['t',0,1,2,3,4,5], [0,0,0,0,0,0,0], [1,0,0,0,0,0,0]]
t_matrix_buy = [['t',0,1,2,3,4,5], [0,0,0,0,0,0,0], [1,0,0,0,0,0,0]]
t_matrix_sell = [['t',0,1,2,3,4,5], [0,0,0,0,0,0,0], [1,0,0,0,0,0,0]]

m_matrix_batteri = [['m',0,1,2,3,4,5], [0,0,0,0,0,0,0], [1,0,0,0,0,0,0]]
m_matrix_H_lagring = [['m',0,1,2,3,4,5], [0,0,0,0,0,0,0], [1,0,0,0,0,0,0]]
m_matrix_buy = [['m',0,1,2,3,4,5], [0,0,0,0,0,0,0], [1,0,0,0,0,0,0]]
m_matrix_sell = [['m',0,1,2,3,4,5], [0,0,0,0,0,0,0], [1,0,0,0,0,0,0]]

vW_matrix_batteri = [['vW',0,1,2,3,4,5], [0,0,0,0,0,0,0], [1,0,0,0,0,0,0]]
vW_matrix_H_lagring = [['vW',0,1,2,3,4,5], [0,0,0,0,0,0,0], [1,0,0,0,0,0,0]]
vW_matrix_buy = [['vW',0,1,2,3,4,5], [0,0,0,0,0,0,0], [1,0,0,0,0,0,0]]
vW_matrix_sell = [['vW',0,1,2,3,4,5], [0,0,0,0,0,0,0], [1,0,0,0,0,0,0]]

vA_matrix_batteri = [['vA',0,1,2,3,4,5], [0,0,0,0,0,0,0], [1,0,0,0,0,0,0]]
vA_matrix_H_lagring = [['vA',0,1,2,3,4,5], [0,0,0,0,0,0,0], [1,0,0,0,0,0,0]]
vA_matrix_buy = [['vA',0,1,2,3,4,5], [0,0,0,0,0,0,0], [1,0,0,0,0,0,0]]
vA_matrix_sell = [['vA',0,1,2,3,4,5], [0,0,0,0,0,0,0], [1,0,0,0,0,0,0]]


for key, value in summary.items():
    input_string = key
    split_string = input_string.split('_')

    if split_string[4] == 'tmvW':
        tmvW_matrix_batteri[int(split_string[6]) + 1][int(split_string[5]) + 1] = value[0]
        tmvW_matrix_H_lagring[int(split_string[6]) + 1][int(split_string[5]) + 1] = value[1]
        tmvW_matrix_buy[int(split_string[6]) + 1][int(split_string[5]) + 1] = value[2]
        tmvW_matrix_sell[int(split_string[6]) + 1][int(split_string[5]) + 1] = value[3]
    elif split_string[4] == 'tmvA':
        tmvA_matrix_batteri[int(split_string[6]) + 1][int(split_string[5]) + 1] = value[0]
        tmvA_matrix_H_lagring[int(split_string[6]) + 1][int(split_string[5]) + 1] = value[1]
        tmvA_matrix_buy[int(split_string[6]) + 1][int(split_string[5]) + 1] = value[2]
        tmvA_matrix_sell[int(split_string[6]) + 1][int(split_string[5]) + 1] = value[3]
    elif split_string[4] == 'tm':
        tm_matrix_batteri[int(split_string[6]) + 1][int(split_string[5]) + 1] = value[0]
        tm_matrix_H_lagring[int(split_string[6]) + 1][int(split_string[5]) + 1] = value[1]
        tm_matrix_buy[int(split_string[6]) + 1][int(split_string[5]) + 1] = value[2]
        tm_matrix_sell[int(split_string[6]) + 1][int(split_string[5]) + 1] = value[3]
    elif split_string[4] == 'mvW':
        mvW_matrix_batteri[int(split_string[6]) + 1][int(split_string[5]) + 1] = value[0]
        mvW_matrix_H_lagring[int(split_string[6]) + 1][int(split_string[5]) + 1] = value[1]
        mvW_matrix_buy[int(split_string[6]) + 1][int(split_string[5]) + 1] = value[2]
        mvW_matrix_sell[int(split_string[6]) + 1][int(split_string[5]) + 1] = value[3]
    elif split_string[4] == 'mvA':
        mvA_matrix_batteri[int(split_string[6]) + 1][int(split_string[5]) + 1] = value[0]
        mvA_matrix_H_lagring[int(split_string[6]) + 1][int(split_string[5]) + 1] = value[1]
        mvA_matrix_buy[int(split_string[6]) + 1][int(split_string[5]) + 1] = value[2]
        mvA_matrix_sell[int(split_string[6]) + 1][int(split_string[5]) + 1] = value[3]
    elif split_string[4] == 'tvW':
        tvW_matrix_batteri[int(split_string[6]) + 1][int(split_string[5]) + 1] = value[0]
        tvW_matrix_H_lagring[int(split_string[6]) + 1][int(split_string[5]) + 1] = value[1]
        tvW_matrix_buy[int(split_string[6]) + 1][int(split_string[5]) + 1] = value[2]
        tvW_matrix_sell[int(split_string[6]) + 1][int(split_string[5]) + 1] = value[3]
    elif split_string[4] == 'tvA':
        tvA_matrix_batteri[int(split_string[6]) + 1][int(split_string[5]) + 1] = value[0]
        tvA_matrix_H_lagring[int(split_string[6]) + 1][int(split_string[5]) + 1] = value[1]
        tvA_matrix_buy[int(split_string[6]) + 1][int(split_string[5]) + 1] = value[2]
        tvA_matrix_sell[int(split_string[6]) + 1][int(split_string[5]) + 1] = value[3]
    elif split_string[4] == 't':
        t_matrix_batteri[int(split_string[6]) + 1][int(split_string[5]) + 1] = value[0]
        t_matrix_H_lagring[int(split_string[6]) + 1][int(split_string[5]) + 1] = value[1]
        t_matrix_buy[int(split_string[6]) + 1][int(split_string[5]) + 1] = value[2]
        t_matrix_sell[int(split_string[6]) + 1][int(split_string[5]) + 1] = value[3]
    elif split_string[4] == 'm':
        m_matrix_batteri[int(split_string[6]) + 1][int(split_string[5]) + 1] = value[0]
        m_matrix_H_lagring[int(split_string[6]) + 1][int(split_string[5]) + 1] = value[1]
        m_matrix_buy[int(split_string[6]) + 1][int(split_string[5]) + 1] = value[2]
        m_matrix_sell[int(split_string[6]) + 1][int(split_string[5]) + 1] = value[3]
    elif split_string[4] == 'vW':
        vW_matrix_batteri[int(split_string[6]) + 1][int(split_string[5]) + 1] = value[0]
        vW_matrix_H_lagring[int(split_string[6]) + 1][int(split_string[5]) + 1] = value[1]
        vW_matrix_buy[int(split_string[6]) + 1][int(split_string[5]) + 1] = value[2]
        vW_matrix_sell[int(split_string[6]) + 1][int(split_string[5]) + 1] = value[3]
    elif split_string[4] == 'vA':
        vA_matrix_batteri[int(split_string[6]) + 1][int(split_string[5]) + 1] = value[0]
        vA_matrix_H_lagring[int(split_string[6]) + 1][int(split_string[5]) + 1] = value[1]
        vA_matrix_buy[int(split_string[6]) + 1][int(split_string[5]) + 1] = value[2]
        vA_matrix_sell[int(split_string[6]) + 1][int(split_string[5]) + 1] = value[3]

matrix_excel_batteri = [tmvW_matrix_batteri, tmvA_matrix_batteri, tm_matrix_batteri, mvW_matrix_batteri, mvA_matrix_batteri, tvW_matrix_batteri, tvA_matrix_batteri, t_matrix_batteri, m_matrix_batteri, vW_matrix_batteri, vA_matrix_batteri]
matrix_excel_H_lagring = [tmvW_matrix_H_lagring, tmvA_matrix_H_lagring, tm_matrix_H_lagring, mvW_matrix_H_lagring, mvA_matrix_H_lagring, tvW_matrix_H_lagring, tvA_matrix_H_lagring, t_matrix_H_lagring, m_matrix_H_lagring, vW_matrix_H_lagring, vA_matrix_H_lagring]
matrix_excel_buy = [tmvW_matrix_buy, tmvA_matrix_buy, tm_matrix_buy, mvW_matrix_buy, mvA_matrix_buy, tvW_matrix_buy, tvA_matrix_buy, t_matrix_buy, m_matrix_buy, vW_matrix_buy, vA_matrix_buy]
matrix_excel_sell = [tmvW_matrix_sell, tmvA_matrix_sell, tm_matrix_sell, mvW_matrix_sell, mvA_matrix_sell, tvW_matrix_sell, tvA_matrix_sell, t_matrix_sell, m_matrix_sell, vW_matrix_sell, vA_matrix_sell]

excel_batteri = pd.DataFrame(matrix_excel_batteri)
excel_H_lagring = pd.DataFrame(matrix_excel_H_lagring)
excel_buy = pd.DataFrame(matrix_excel_buy)
excel_sell = pd.DataFrame(matrix_excel_sell)

# Write DataFrames to Excel
with pd.ExcelWriter('test123.xlsx') as writer:
    excel_batteri.to_excel(writer, sheet_name='batteri', index=False, header=False)
    excel_H_lagring.to_excel(writer, sheet_name='H_lagring', index=False, header=False)
    excel_buy.to_excel(writer, sheet_name='buy_electricity', index=False, header=False)
    excel_sell.to_excel(writer, sheet_name='sell_electricity', index=False, header=False)
#-------------------------------------------------------------------------------------------


intressanta = ['result_sum_power_tmvW_1_1', 'result_sum_power_tmvA_1_1', 'result_sum_power_tm_1_1']



for plotta in intressanta:
    batteri = results[plotta][0]
    H_lagring = results[plotta][1]
    buy_electricity = results[plotta][2]
    sell_electricity = results[plotta][3]

    plt.subplot(4, 1, 1)  # (rows, columns, index)
    plt.plot(range(len(sum_power) + 1), batteri)  # Plot your data for the first subplot
    plt.title('Laddning Batteri '+plotta)
    plt.xlabel('Tid (h)')
    plt.ylabel('Energi (kWh)')

    # Subplot 2
    plt.subplot(4, 1, 2)
    plt.plot(range(len(sum_power) + 1), H_lagring)  # Plot your data for the second subplot
    plt.title('Laddning Vätgas '+plotta)
    plt.xlabel('Tid (h)')
    plt.ylabel('Energi (kWh)')

    # Subplot 3
    plt.subplot(4, 1, 3)
    plt.plot(range(len(sum_power) + 1), buy_electricity)  # Plot your data for the third subplot
    plt.title('Köpt energi '+plotta)
    plt.xlabel('Tid (h)')
    plt.ylabel('Energi (kWh)')

    plt.subplot(4, 1, 4)
    plt.plot(range(len(sum_power) + 1), sell_electricity)  # Plot your data for the third subplot
    plt.title('Såld energi '+plotta)
    plt.xlabel('Tid (h)')
    plt.ylabel('Energi (kWh)')

    plt.tight_layout()  # Adjust layout to prevent overlapping
    plt.show()