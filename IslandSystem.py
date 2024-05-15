import SpecsMalsta as Specs
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import openpyxl
from datetime import datetime

def data_extender(data_list, years):
    data_extended = []
    data = data_list
    data_extended.extend(data)
    for i in range(years-1):
        data_extended.extend(data)
    return data_extended

def system(sum_power,x, y):
    batteri_spec, H_lagring_spec = Specs.specifikationer(x, y)

    batteri = np.zeros(len(sum_power) + 1, dtype=float).tolist()
    batteri[0] = batteri_spec[0] * 1  # fulladdat batteri, kan multipliceras till fler batterier

    H_lagring = np.zeros(len(sum_power) + 1, dtype=float).tolist()
    # H_lagring[0] = H_lagring_spec[0]*1 # full tank vätgas, kan multipliceras för fler tankar
    H_lagring[0] = 0  # tom tank i början

    buy_electricity = np.zeros(len(sum_power) + 1, dtype=float).tolist()
    sell_electricity = np.zeros(len(sum_power) + 1, dtype=float).tolist()

    for i in range(len(sum_power)):  # Gå över alla punkter i sum_power
        # ------------ Om ingen energi stoppas in eller tas ur vätgaslagring eller batteri så behåller de sin laddning
        batteri[i + 1] = batteri[i]
        H_lagring[i + 1] = H_lagring[i]

        # ------------ Block för när energin producerad är mindre än den som konsumeras
        if sum_power[i] < 0 and batteri[i] == 0 and H_lagring == 0:  # Om produktion är mindre än konsumtion och batteriet samt vätgaslagring är tomt
            buy_electricity[i] = -sum_power[i]  # Köper el, - för att vi vill ha positivt värde på W vi köper
        elif sum_power[i] < 0 and batteri[i] == 0:  # Om produktion är mindre än konsumtion och batteriet är slut men det finns laddning kvar i vätgaslagringen
            H_lagring_step = H_lagring[i] + sum_power[i]  # Vill se hur mycket energi vi kan ta ut, kikar med ett steg framåt. + pga negativ sum_power
            if H_lagring_step < 0:  # Om vi tar ut så mycket energi vi vill ha kommer det resultera i en tom vätgastank
                H_lagring[i + 1] = 0  # Tom vätgastank
                buy_electricity[i] = -H_lagring_step  # Den el vi måste köpa. + pga negativ sum_power
            else:
                H_lagring[i + 1] = H_lagring[i] + sum_power[i]  # Om inte, ta ut energi som vanligt från vätgaslagret. + pga negativ sum_power
        elif sum_power[i] < 0:  # Om produktion är mindre än konsumtion men det finns energi i batteriet och vätgaslagret
            batteri_step = batteri[i] + sum_power[i]  # Se ett steg framåt om laddningen kommer ta slut i batteriet
            H_lagring_step = batteri[i] + H_lagring[i] + sum_power[i]  # Se ett steg framåt om laddning för både batteri och vätgas kommer att ta slut
            if batteri_step == 0:  # Om batteriet tar exakt slut behövs ingen mer energi tas från annan plats
                batteri[i + 1] = 0  # sätt laddning till noll
            elif batteri_step < 0:  # Om batteriet kommer att ta mer än slut
                batteri[i + 1] = 0  # sätt laddning till noll
                if H_lagring_step < 0:  # Om vätgaslagret även kommer att ta slut
                    H_lagring[i + 1] = 0  # Sätt laddning till noll
                    buy_electricity[i] = -H_lagring_step  # köp resterande energi. - På grund av H_lagring_step har negativt värde
                else:
                    H_lagring[i + 1] = H_lagring[i] + sum_power[i]  # Om vätgaslagret inte kommer att ta slut, eller precis ta slut, beräkna hur mycket som är kvar
            else:
                batteri[i + 1] = batteri[i] + sum_power[i]  # Om det kommer att finnas laddning kvar i batteriet.+ för att sum_power är negativt

        if sum_power[i] > 0 and batteri[i] == batteri_spec[0] and H_lagring[i] == H_lagring_spec[0]:  # Om produktion är större än konsumtion och både batteri är vätgaslagring är fulla
            sell_electricity[i] = sum_power[i]  # Sälj elektriciteten
        elif sum_power[i] > 0 and batteri[i] == batteri_spec[0]:  # Om produktion är större än konsumtion och endast batteriet är fulladdat
            H_lagring_step = H_lagring[i] + Specs.lagring_H(sum_power[i])  # Kolla ett steg framåt för att se om energin överstiger maxladdning för vätgastank. Skickar sum_power till lagring_h för beräkning av förluster
            if H_lagring_step > H_lagring_spec[0] or H_lagring_step == H_lagring_spec[0]:  # Om energin kommer överfylla tanken, toppa upp tanken och sälj resterande energi
                H_lagring[i + 1] = H_lagring_spec[0]  # fyller tanken
                sell_electricity[i] = H_lagring[i] + sum_power[i] - H_lagring_spec[0]  # säljer resterande energi
            else:
                H_lagring[i + 1] = H_lagring[i] + Specs.lagring_H(sum_power[i])  # Annars, fyll tank med energi
        elif sum_power[i] > 0:  # Om produktion är större än konsumtion men batteri och vätgas är inte fulladdade
            batteri_step = batteri[i] + Specs.batteri(sum_power[i])  # Se ett steg framåt
            if batteri_step > batteri_spec[0]:  # Om batteriet blir fulladdat
                batteri[i + 1] = batteri_spec[0]  # ladda batteriet
                H_lagring_step = H_lagring[i] + Specs.lagring_H(sum_power[i])
                if H_lagring_step > H_lagring_spec[0]:
                    H_lagring[i + 1] = H_lagring_spec[0]
                    sell_electricity[i] = H_lagring[i] + sum_power[i] - H_lagring_spec[0]  # säljer resterande energi
                else:
                    H_lagring[i + 1] = H_lagring[i] + batteri_step - batteri_spec[0]  # Resten går in i vätelagring
            elif batteri_step < batteri_spec[0] or batteri_step == batteri_spec[0]:  # Om batteriet inte blir fullt eller om det blir exakt full
                batteri[i + 1] = batteri[i] + Specs.batteri(sum_power[i])  # Ladda batteriet

        if sum_power[i] == 0:  # Perfekt, gör inget
            pass

    return batteri, H_lagring, buy_electricity, sell_electricity