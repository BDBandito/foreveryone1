import math
import matplotlib.pyplot as plt
import numpy as np
from decimal import Decimal, getcontext
import csv
getcontext().prec = 40 # Anger max antal decimaler för Decimal

# transienter som används
transients = ['e', 'N_u238', 'alfa', 'andel_PuvsU', 'Sigma_a235', 'Sigma_f235', 'Sigma_a238', 'Sigma_a239', 'Sigma_f239',
              'Sigma_abransle', 'p', 'C', 'Nt_forbrukade', 'Nt_prod', 'Nt_fissila', 'Nt_235', 'Nt_239', 'Nt_238', 'ny',
              'fi_bransle', 'gamma_I', 'gamma_Xe', 'N_I', 'm_I', 'N_Xe', 'm_Xe', 'V_gift', 'Sigma_agift', 'eta', 'f', 'k']

# Konstanter för UO2 eller UN. Av-/kommenteras vid val av bränsle
# För UN:
rho_U = Decimal(14.3)  # [g/cm3] bränslets densitet
T_u = Decimal(500 + 273.15)  # [K] Bränsletemperatur
atomvikt_U235 = 235 + 15 # atomvikt för bränslet med U235
atomvikt_U238 = 238 + 15 # atomvikt för bränslet med U238
R_u = Decimal(0.43)  # [cm] stavradie
#####
# För UO2
#rho_U = Decimal(10.5)  # [g/cm3] bränslets densitet
#T_u = Decimal(1000 + 273.15)  # [K] Bränsletemperatur
#atomvikt_U235 = 235 + 16*2 # atomvikt för bränslet med U235
#atomvikt_U238 = 238 + 16*2 # atomvikt för bränslet med U238
#R_u = Decimal(0.5)  # [cm] stavradie
#####


# Konstanter
e_start = Decimal(0.020791)  # anrikningsvärde vid start
epsilon = Decimal(1.06)  # Snabba fissionsfaktorn
P_leakage = Decimal(0.97)  # Icke-läckagefaktorn
P_s = Decimal(0.96)  # icke-läckagefaktor för snabba neutroner - taget från första bästa källa på nätet
P = Decimal(3300 * 10 ** 6)  # [W] Termisk reaktoreffekt
m_UO2 = Decimal(122000)  # [kg]vikt för uranoxidet
u = Decimal(1.66 * 10 ** (-24))  # [g] universala massenheten
E_fission = Decimal(1.602 * 10 ** (-19) * 200 * 10 ** 6)  # [J] utlöst energi per fission
ny_pu239 = Decimal(2.9)  # genomsnittligt antal för hur många neutroner som bildas per fission för Pu239
ny_u235 = Decimal(2.42)  # genomsnittligt antal för hur många neutroner som bildas per fission för U235
sigma_f235 = Decimal(582.2 * 10 ** -24)  # [cm^2] mikroskopiskt tvärsnitt för fission U235
sigma_c235 = Decimal(93.6 * 10 ** -24)  # [cm^2] mikroskopiskt tvärsnitt för infångning U235
sigma_c238 = Decimal(2.7 * 10 ** -24)  # [cm^2] mikroskopiskt tvärsnitt för infångning U238
sigma_a238 = Decimal(2.7 * 10 ** -24) # [cm^2] mikroskopiskt tvärsnitt för absorption U238
sigma_f238 = Decimal(16.8 * 10 ** (-6) * 10 ** (-24)) # [cm^2]
N0_u235 = (e_start * rho_U) / ((Decimal(atomvikt_U235) * u))  # [kärnor/cm3] kärntäthet U235 vid start
N0tot_u235 = N0_u235 * (m_UO2 / (rho_U * Decimal(10 ** -3)))  # antal kärnor U235 totalt i härden vid start
N0_u238 = ((1 - e_start) * rho_U) / ((atomvikt_U235) * u)  # [kärnor/cm3] kärntäthet U238 vid start, KSU s102
N0tot_u238 = N0_u238 * (m_UO2 / (rho_U * Decimal(10 ** -3)))  # antal kärnor U238 totalt i härden vid start
e_start = N0tot_u235/(N0tot_u238+N0tot_u235)
N0_pu239 = 0  # antal Pu239 kärnor vid start
sigma_a235 = sigma_f235 + sigma_c235  # [cm^2] Mikroskopiskt tvärsnitt för absorption för U235
sigma_f239 = Decimal(744.4 * 10 ** -24)  # [cm^2] Mikroskopiskt fissionstvärsnitt för Pu239
sigma_c239 = Decimal(268.8 * 10 ** -24)  # [cm^2] Mikroskopiskt infångningstvärsnitt för Pu239
sigma_a239 = sigma_f239 + sigma_c239 # [cm^2] Mikroskopiskt absorptionstvärsnitt för Pu239
rho_Xe = 5.9  # [g/cm3] densitet för Xe135
FR = P / (Decimal(200 * 1.6021 * (10 ** -13)))  #[1/cm3/s] fissionsrat för bränslet

V_u = Decimal(1)  # [geometriskt förhållande] Volym bränsle
V_m = Decimal(3)  # [geometriskt förhållande] Volym moderator
xi = Decimal(1.26)  # log. medelenergiförlusten för moderatorn
Sigma_s = Decimal(1)  # [cm-1] makroskopiskt spridningstvärsnitt i moderatormaterialet


sigma_Iabs_300K = 3 + Decimal(39.6 / (math.sqrt(R_u * rho_U))) # [barn] Resonansintegral vid 300K
B_1 = Decimal(6.1 * 10 ** (-3)) + (Decimal(0.94 * 10 ** (-2)) / (R_u * rho_U))
sigma_IabsT = sigma_Iabs_300K * (1 + B_1 * Decimal((math.sqrt(T_u) - math.sqrt(300))))  # [barn] Bränslets resonansintegral, används vid beräkning av resonanspassagefaktorn

insattning = Decimal(0)  # [%] insättning av styrstaven
V_bransle = Decimal(1)  # typiskt volymförhållande

fi_mod = Decimal(1) # neutronflöde moderator
Sigma_amod = Decimal(0.02)  # taget från KSU (I think)
V_mod = Decimal(3)  # typiskt volymförhållande

fi_styrstav = Decimal(1) # neutronflöde moderator
Sigma_astyrstav = Decimal(0.15)  # [cm-1]
styrstav_bredd = Decimal(270 * 10 ** -3)  # [m]
styrstav_tjocklek = Decimal(7 * 10 ** -3)  # [m]
styrstav_langd = Decimal(3700 * 10 ** -3)  # [m]
N_styrstav = Decimal(161)  # antal styrstavar
V_styrstav = styrstav_bredd * styrstav_tjocklek * styrstav_langd * insattning * N_styrstav  # Volym styrstav

fi_gift = Decimal(1) # Neutronflöde reaktorgift
sigma_agift = Decimal(2.65 * 10 ** 6)  # [barn] mikroskopiskt tvärsnitt för Xe135, taget från KSU
gamma_TeU = Decimal(0.0321618)  # fissionsproduktfördelning för Tellurium från U235
gamma_IU = Decimal(0.0292737)  # fissionsproduktfördelning för jod från U235
T_I2 = Decimal(6.57 * 60 * 60)  # [s] halveringstid för jod
lambda_I = Decimal(math.log(2)) / T_I2  # sönderfallskonstanten för I135

gamma_TePu = Decimal(0.0219297)  # fission yield Tellurium från Pu239
gamma_IPu = Decimal(0.0428742)  # fission yield Iodine från Pu239
gamma_XeU = Decimal(0.00178122)  # fission yield Xenon från U235
gamma_XePu = Decimal(0.007522799)  # fission yield Xenon från Pu239
N_A = Decimal(6.022 * 10 ** 23)  # Avogadros konstant
M_I = Decimal(134.91006)  # [u] relativ atommassa
T_Xe2 = Decimal(9.14 * 60 * 60)  # [s] halveringstid Xe135
lambda_Xe = Decimal(math.log(2)) / T_Xe2  # [s] sönderfallskonstant för Xe135
sigma_aXe = Decimal(2.65 * 10 ** 6 * 10 ** -24)  # [cm^2] mikroskopiskt absorptionstvärsnitt för Xe135
M_Xe = Decimal(134.9072075)  # [u]Atomic mass Xe135

fi_an = Decimal(1) # neutronflöde annat
Sigma_aan = Decimal(1) # [cm^-1] Makroskopiskt absorptionstvärsnitt annat
V_an = Decimal(1) # volym annat


# Listor för att lagra transienter
e = []
N_u238 = []
alfa = []
andel_PuvsU = []
Sigma_a235 = []
Sigma_a238 = []
Sigma_f235 = []
Sigma_a239 = []
Sigma_f239 = []
Sigma_abransle = []
p = []
C = []
Nt_forbrukade = []
Nt_prod = []
Nt_fissila = []
Nt_235 = []
Nt_239 = []
Nt_238 = []
ny = []

fi_bransle = []
gamma_I = []
gamma_Xe = []
N_I = []
m_I = []
N_Xe = []
m_Xe = []
V_gift = []
Sigma_agift = []
eta = []
f = []
k = []

e_out = []
N_u238_out = []
alfa_out = []
andel_PuvsU_out = []
Sigma_a235_out = []
Sigma_a238_out = []
Sigma_f235_out = []
Sigma_a239_out = []
Sigma_f239_out = []
Sigma_abransle_out = []
p_out = []
C_out = []
Nt_forbrukade_out = []
Nt_prod_out = []
Nt_fissila_out = []
Nt_235_out = []
Nt_239_out = []
Nt_238_out = []
ny_out = []

fi_bransle_out = []
gamma_I_out = []
gamma_Xe_out = []
N_I_out = []
m_I_out = []
N_Xe_out = []
m_Xe_out = []
V_gift_out = []
Sigma_agift_out = []
eta_out = []
f_out = []
k_out = []

# Sätter startvärden
e.insert(0, e_start) # Anrikning
N_u238.insert(0, N0_u238) # [kärnor/cm^3]Kärntäthet U238
alfa.insert(0, Decimal(0.5892255555)) # void
alfa.insert(1, Decimal(0.5892255555)) # void
Sigma_a235.insert(0, ((e_start * rho_U) / (235 * u)) * sigma_a235) # [cm^-1] makroskopiska tvärsnittet för absorption för uran235
Sigma_a238.insert(0, (((1 - e_start) * rho_U) / (238 * u)) * sigma_a238)  # [cm^-1] makroskopiska tvärsnittet för absorption för uran238
Sigma_f235.insert(0, (((e_start) * rho_U) / (235 * u)) * sigma_f235)  # [cm^-1] makroskopiska tvärsnittet för fission för uran235
Sigma_a239.insert(0, 0) # [cm^-1] makroskopiska tvärsnittet för absorption för Pu239
Sigma_f239.insert(0, 0) # [cm^-1] makroskopiska tvärsnittet för fission för Pu239
Sigma_abransle.insert(0, Sigma_a235[0] + Sigma_a238[0] + Sigma_a239[0])  # [cm^-1] Makroskopiska tvärsnittet för bränslet, ändras över tiden
p.insert(0, Decimal(math.exp(-(1 - e_start) * ((N0_u238 * sigma_IabsT * Decimal((10 ** -24)) * V_u) / (xi * Sigma_s * V_m * (1 - alfa[0])))))) # resonanspassagefaktorn
C.insert(0, Decimal((Sigma_a238[0] / Sigma_a235[0]) + ((Sigma_f235[0] / (Sigma_a235[0] ) * epsilon * ny_u235 * (1 - p[0]))))) # konverteringskvoten

DeltaNt_forbrukade = 0 # Antal förbrukade kärnor detta steg. Inga till en början
Nt_forbrukade.insert(0, 0)  # antal förbrukade fissila kärnor totalt. Inga till en början

DeltaNt_prod = 0 # Antal producerade kärnor detta steg. Inga till en början
Nt_prod.insert(0, 0) # Antal producerade kärnor. Inga till en början
Nt_fissila.insert(0, N0tot_u235) # Antal fissila kärnor
Nt_235.insert(0, Decimal(N0tot_u235))  # Antal kärnor U235
Nt_239.insert(0, 0)  # Antal kärnor Pu239
Nt_238.insert(0, N0tot_u238)  # Antal kärnor U238
andel_PuvsU.insert(0, 0) # Andel PuvsU, noll till en början för att det inte finns Pu239
N_u238.insert(0, ((1 - e[0]) * rho_U) / ((atomvikt_U238) * u)) # [kärnor/cm^3] Kärntäthet U238


ny.insert(0, ny_u235) # Antal nya neutroner per klyvning
fi_bransle.insert(0, Sigma_abransle[0] / FR)  # neutronflöde kopplat till absorption för U235
gamma_I.insert(0, (gamma_TeU + gamma_IU))  # total fission yield för Jod, när vi vet att det finns Pu239 i reaktorn
gamma_Xe.insert(0, gamma_XeU)  # total fission yield för Xenon, när vi vet att det finns Pu239 i reaktorn
N_I.insert(0, (gamma_I[0] * FR) / lambda_I)  # antalet kärnor Jod i härden
m_I.insert(0, (N_I[0] / N_A) * M_I)  # [g] massa jod i härden
N_Xe.insert(0, ((gamma_Xe[0] + gamma_I[0]) * FR) / (gamma_Xe[0] + sigma_aXe * fi_bransle[0]))  # antal kärnor Xe135 i härden - denna tagen från slide (är nog denna)
m_Xe.insert(0, (N_Xe[0] / N_A) * M_Xe)  # [g] massa Xe135 i härden
V_gift.insert(0, N_Xe[0] / ((N0tot_u235 + N0tot_u238) * 2 * 16))  # Volym reaktorgift
Sigma_agift.insert(0, N_Xe[0] * sigma_aXe)  # Makroskopiska tvärsnittet för reaktorgift

eta.insert(0, ny[0] * (e[0] * sigma_f235) / (e[0] * (sigma_a235) + (1 - e[0]) * sigma_c238)) # Termiska fissionsfaktorn
f.insert(0, (Sigma_abransle[0] * V_bransle) / (Sigma_abransle[0] * V_bransle + Sigma_amod * V_mod * (1 - alfa[0]) + Sigma_astyrstav * V_styrstav + Sigma_agift[0] * V_gift[0] + Decimal(0.15))) # Termiska utnyttjandefaktorn
k.insert(0, eta[0] * epsilon * p[0] * f[0] * P_leakage) # Multiplikationsfaktorn

stop_time = [] # tom lista, denna fylls när det är dags att stoppa simulationen

def next_step(past, present):
    Sigma_a235.insert(present, ((e[past] * (1 - andel_PuvsU[past]) * rho_U) / (235 * u)) * sigma_a235) # [cm^-1] makroskopiskt tvärsnitt för absorption för U235
    Sigma_a238.insert(present, (((1 - e[past]) * rho_U) / (238 * u)) * sigma_a238)  # [cm^-1] makroskopiska tvärsnittet för absorption för uran238
    Sigma_f235.insert(present, (((e[past] * (1 - andel_PuvsU[past])) * rho_U) / (235 * u)) * sigma_f235)  # [cm^-1] makroskopiska tvärsnittet för fission för uran235
    Sigma_a239.insert(present, ((e[past] * andel_PuvsU[past] * rho_U) / (239 * u)) * sigma_a239) # [cm^-1] makroskopiska tvärsnittet för absorption för uran235
    Sigma_f239.insert(present, ((e[past] * andel_PuvsU[past] * rho_U) / (239 * u)) * sigma_f239) # [cm^-1] makroskopiska tvärsnittet för fission för Pu239
    Sigma_abransle.insert(present, Sigma_a235[past] + Sigma_a238[past] + Sigma_a239[past])  # [cm^-1]Makroskopiska tvärsnittet för bränslet, ändras över tiden
    C.insert(present, Decimal((Sigma_a238[past] / (Sigma_a235[past] * (1 - andel_PuvsU[past]) + Sigma_a239[past] * andel_PuvsU[past])) + ((Sigma_f235[past] * (1 - andel_PuvsU[past]) + Sigma_f239[past] * andel_PuvsU[past]) / (Sigma_a235[past] * (1 - andel_PuvsU[past]) + Sigma_a239[past] * andel_PuvsU[past])) * epsilon * ny_u235 * (1 - p[past]))) # konverteringskvot

    DeltaNt_forbrukade = ((P * timestep) / E_fission) * Decimal(0.94) * ((sigma_a235 / sigma_f235) * (1 - andel_PuvsU[past]) + (sigma_a239 / sigma_f239) * andel_PuvsU[past])  # Gammal formel
    Nt_forbrukade.insert(present, DeltaNt_forbrukade + Nt_forbrukade[past])  # antal förbrukade fissila kärnor totalt

    DeltaNt_prod = DeltaNt_forbrukade * C[past] # Antal producerade kärnor detta steg
    Nt_prod.insert(present, Decimal(DeltaNt_prod + Nt_prod[past])) # Antal producerade kärnor totalt
    Nt_fissila.insert(present, Decimal(Nt_239[past] + Nt_235[past])) # Totalt antal fissila kärnor

    e.insert(present, (Nt_235[past] + Nt_239[past]) / (Nt_235[past] + Nt_239[past] + Nt_238[past])) # anrikning
    kvot_Sigma_239_235 = (sigma_f239 * Nt_239[past]) / ((sigma_f235 * Nt_235[past]) + (sigma_f239 * Nt_239[past])) # Kvot mellan makroskopiska fissionstvärsnitt
    Nt_235.insert(present, Decimal(N0tot_u235 - Decimal(0.94) * Nt_forbrukade[past] * (1 - kvot_Sigma_239_235)))  # Antal kärnor U235
    Nt_239.insert(present, Decimal(Nt_prod[past] - Decimal(0.94) * Nt_forbrukade[past] * kvot_Sigma_239_235))  # Antal kärnor Pu239
    Nt_238.insert(present, N0tot_u238 - Decimal(0.06)*Nt_forbrukade[past])  # Antal kärnor U238
    andel_PuvsU.insert(present, Nt_239[past] / (Nt_239[past] + Nt_235[past])) # Andel PuvsU
    N_u238.insert(present, ((1 - e[past]) * rho_U) / (atomvikt_U238 * u)) # [kärnor/cm^3] kärntäthet U238


    ny.insert(present, (1 - andel_PuvsU[past]) * ny_u235 + (andel_PuvsU[past]) * ny_pu239) # Antal frisläppta neutroner per klyvning
    fi_bransle.insert(present, Sigma_abransle[past] / FR)  # neutronflöde kopplat till absorption för U235
    gamma_I.insert(present, (gamma_TeU + gamma_IU) * (1 - andel_PuvsU[1]) + (gamma_TePu + gamma_IPu) * andel_PuvsU[past])  # total fission yield för Jod, när vi vet att det finns Pu239 i reaktorn
    gamma_Xe.insert(present, gamma_XeU * (1 - andel_PuvsU[past]) + gamma_XePu * andel_PuvsU[past])  # total fission yield för Xenon, när vi vet att det finns Pu239 i reaktorn
    N_I.insert(present, (gamma_I[past] * FR) / lambda_I)  # antalet kärnor Jod i härden
    m_I.insert(present, (N_I[past] / N_A) * M_I)  # [g] massa jod i härden
    N_Xe.insert(present, ((gamma_Xe[past] + gamma_I[past]) * FR) / (gamma_Xe[past] + sigma_aXe * fi_bransle[past]))  # antal kärnor Xe135 i härden - denna tagen från slide (är nog denna)
    m_Xe.insert(present, (N_Xe[past] / N_A) * M_Xe)  # [g] massa Xe135 i härden
    V_gift.insert(present, N_Xe[past] / ((N0tot_u235 + N0tot_u238) * 2 * 16)) # volym reaktorgift
    Sigma_agift.insert(present, N_Xe[past] * sigma_aXe)  # Makroskopiska tvärsnittet för reaktorgift, vilket i detta fall endast är Xe-135

    eta.insert(present, ny[past] * (e[past] * (sigma_f235 * (1 - andel_PuvsU[past]) + sigma_f239 * andel_PuvsU[past]) / (e[past] * (sigma_a235 * (1 - andel_PuvsU[past]) + sigma_a239 * andel_PuvsU[past]) + (1 - e[past]) * sigma_c238))) # termiska fissionsfaktorn
    p.insert(present, Decimal(math.exp(-(1 - e[past]) * ((N_u238[past] * sigma_IabsT * Decimal((10 ** -24)) * V_u) / (xi * Sigma_s * V_m * (1 - alfa[past])))))) # resonanspassagefaktorn
    f.insert(present, (Sigma_abransle[past] * V_bransle) / (Sigma_abransle[past] * V_bransle + Sigma_amod * V_mod * (1 - alfa[past]) + Sigma_astyrstav * V_styrstav + Sigma_agift[past] * V_gift[past] + Decimal(0.15))) # termiska utnyttjandefaktorn
    k.insert(present, eta[1] * epsilon * p[1] * f[1] * P_leakage) # multiplikationsfaktorn

    # If-satserna och whole-looperna nedanför används för att tillsammans reglera voiden.
    # Om k<1 så sänks voiden. Om k<1 höjs voiden. Om k<0.1 så avslutas simulationen.
    # Om voiden är 0.15 så sänks den inte mer. Detta för att simulera max flöde för HC-pump.
    if k[present] < 1: # Om true - Voiden sänks
        print('counting -alfa')
        if k[present] < 0.1: # Om True - stoppas simulationen
            print('Stopp pga k<0.1')
            stop_time.append('stop')


        alfa_alfa = alfa[present]
        k_alfa = k[present]
        while k_alfa < 1: # Loopen sänker voiden tills att k<1

            alfa_alfa = alfa_alfa - Decimal(0.00001)
            if alfa_alfa < 0.15: # Sänker inte voiden mer än 0.15 pga maximum HC-flöde
                alfa.insert(present, alfa_alfa) # Ger nya voiden
                break

            p_alfa = Decimal(math.exp(-(1 - e[present]) * ((N_u238[present] * sigma_IabsT * Decimal((10 ** -24)) * V_u) / (xi * Sigma_s * V_m * (1 - alfa_alfa)))))
            k_alfa = eta[present] * epsilon * p_alfa * f[present] * P_leakage
            if k_alfa > 1:
                alfa.insert(present, alfa_alfa) # Ger nya voiden


    elif k[present] > 1.0001 and alfa[present] < 0.8: # Höjer voiden när k > 1.0001. Väldigt användbart vid start av simulationen om k inte är exakt 1
        alfa_alfa = alfa[present]
        k_alfa = k[present]
        while k_alfa > 1.0001:

            alfa_alfa = alfa_alfa + Decimal(0.00001)
            p_alfa = Decimal(math.exp(-(1 - e[present]) * ((N_u238[present] * sigma_IabsT * Decimal((10 ** -24)) * V_u) / (xi * Sigma_s * V_m * (1 - alfa_alfa)))))
            k_alfa = eta[present] * epsilon * p_alfa * f[present] * P_leakage
            if k_alfa < 1.0001:
                alfa.insert(present, alfa_alfa) # Ger nya voiden

    else:
        alfa.insert(present, alfa[past]) # Om voiden inte behöver ändras så behålls nuvarande void




timestep = 60 # Hur många sekunder per steg, alltså 1 minut passerar per steg
switch = True
t = 0 # startvärde för tid
list_iterate = 1 # startvärde för list_iterate
t_max = 100000 # Högst antal steg tills simuleringen slutar

# While loopen kör simuleringen tills antingen
while t < t_max and len(stop_time) == 0:

    # If-satsen byter plats på vilken plats som listan ska läsas ifrån.
    # Den senaste platsen som är skriven på används för att beräkna det nya steget.
    # Det tidigare steget skrivs sedan över med ny data.
    if switch:
       next_step(0,1)
    else:
       next_step(1,0)
    t += 1
    switch = not switch

    # If-satsen ser till att endast det första steg och var tusende tidssteg sparas.
    # Då vi är endast är intresserade av beteendet över en hel bränslecykel
    # behöver vi inte spara alla punkter. Genom att spara ett fåtal punkter
    # behåller vi möjligheten att plotta släta grafer utan att ställa för stora
    # krav på hur mycket data som datorn behöver hålla i.
    if list_iterate % 1000 == 0:
        write = True
    elif list_iterate == 1:
        write = True
    else:
        write = False
    if write:
        print(t) # Skriver ut vilket tidssteg simuleringen är på, kommentera för att öka hastighet för simuleringen
        e_out.append(e[1])
        e_out.append(e[0])

        N_u238_out.append(N_u238[1])
        N_u238_out.append(N_u238[0])

        alfa_out.append(alfa[1])
        alfa_out.append(alfa[0])

        andel_PuvsU_out.append(andel_PuvsU[1])
        andel_PuvsU_out.append(andel_PuvsU[0])

        Sigma_a235_out.append(Sigma_a235[1])
        Sigma_a235_out.append(Sigma_a235[0])

        Sigma_a238_out.append(Sigma_a238[1])
        Sigma_a238_out.append(Sigma_a238[0])

        Sigma_f235_out.append(Sigma_f235[1])
        Sigma_f235_out.append(Sigma_f235[0])

        Sigma_a239_out.append(Sigma_a239[1])
        Sigma_a239_out.append(Sigma_a239[0])

        Sigma_f239_out.append(Sigma_f239[1])
        Sigma_f239_out.append(Sigma_f239[0])

        Sigma_abransle_out.append(Sigma_abransle[1])
        Sigma_abransle_out.append(Sigma_abransle[0])

        p_out.append(p[1])
        p_out.append(p[0])

        C_out.append(C[1])
        C_out.append(C[0])

        Nt_forbrukade_out.append(Nt_forbrukade[1])
        Nt_forbrukade_out.append(Nt_forbrukade[0])

        Nt_prod_out.append(Nt_prod[1])
        Nt_prod_out.append(Nt_prod[0])

        Nt_fissila_out.append(Nt_fissila[1])
        Nt_fissila_out.append(Nt_fissila[0])

        Nt_235_out.append(Nt_235[1])
        Nt_235_out.append(Nt_235[0])

        Nt_239_out.append(Nt_239[1])
        Nt_239_out.append(Nt_239[0])

        Nt_238_out.append(Nt_238[1])
        Nt_238_out.append(Nt_238[0])

        ny_out.append(ny[1])
        ny_out.append(ny[0])

        fi_bransle_out.append(fi_bransle[1])
        fi_bransle_out.append(fi_bransle[0])

        gamma_I_out.append(gamma_I[1])
        gamma_I_out.append(gamma_I[0])

        gamma_Xe_out.append(gamma_Xe[1])
        gamma_Xe_out.append(gamma_Xe[0])

        N_I_out.append(N_I[1])
        N_I_out.append(N_I[0])

        m_I_out.append(m_I[1])
        m_I_out.append(m_I[0])

        N_Xe_out.append(N_Xe[1])
        N_Xe_out.append(N_Xe[0])

        m_Xe_out.append(m_Xe[1])
        m_Xe_out.append(m_Xe[0])

        V_gift_out.append(V_gift[1])
        V_gift_out.append(V_gift[0])

        Sigma_agift_out.append(Sigma_agift[1])
        Sigma_agift_out.append(Sigma_agift[0])

        eta_out.append(eta[1])
        eta_out.append(eta[0])

        f_out.append(f[1])
        f_out.append(f[0])

        k_out.append(k[1])
        k_out.append(k[0])

    list_iterate += 1

# Plottar all data för att få en överblick över hur simuleringen gick. Vid kortare tidssteg -> fler data
# Kan därför vara nödvändigt att mutea denna del av koden. Detta pga för lite ram för att hålla all data ska plottas
x_plot = np.linspace(0, t, len(k_out))

plt.subplot(4,3,1)
plt.plot(x_plot, k_out)
plt.plot(x_plot, alfa_out)

plt.subplot(4,3,2)
plt.plot(x_plot, alfa_out)

plt.subplot(4,3,3)
plt.plot(x_plot, Nt_235_out)
plt.plot(x_plot, Nt_239_out)

plt.subplot(4,3,4)
plt.plot(x_plot, Nt_239_out)

plt.subplot(4,3,5)
plt.plot(x_plot, Nt_238_out)

plt.subplot(4,3,6)
plt.plot(x_plot, Nt_fissila_out)

plt.subplot(4,3,7)
plt.plot(x_plot, Nt_forbrukade_out)

plt.subplot(4,3,8)
plt.plot(x_plot, e_out)

plt.subplot(4,3,9)
plt.plot(x_plot, C_out)

plt.subplot(4,3,10)
plt.plot(x_plot, p_out)

plt.subplot(4,3,11)
plt.plot(x_plot, eta_out)

plt.subplot(4,3,12)
plt.plot(x_plot, ny_out)

plt.show()
print(list_iterate) # skriver ut vilken minut som simuleringen avslutar på


#Resterande delen av koden sparar all data i en CSV-fil

data = zip(e_out, N_u238_out, alfa_out, andel_PuvsU_out, Sigma_a235_out, Sigma_f235_out, Sigma_a238_out, Sigma_a239_out, Sigma_f239_out,Sigma_abransle_out, p_out, C_out, Nt_forbrukade_out, Nt_prod_out, Nt_fissila_out, Nt_235_out, Nt_239_out, Nt_238_out, ny_out, fi_bransle_out,gamma_I_out, gamma_Xe_out, N_I_out, m_I_out, N_Xe_out, m_Xe_out, V_gift_out, Sigma_agift_out, eta_out, f_out, k_out)

file_name = 'my_data2.csv' # namn till filen

# skriver över all kod till filen
with open(file_name, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['e', 'N_u238', 'alfa', 'andel_PuvsU', 'Sigma_a235', 'Sigma_f235', 'Sigma_a238', 'Sigma_a239','Sigma_f239','Sigma_abransle', 'p', 'C', 'Nt_forbrukade', 'Nt_prod', 'Nt_fissila', 'Nt_235', 'Nt_239', 'Nt_238', 'ny','fi_bransle','gamma_I', 'gamma_Xe', 'N_I', 'm_I', 'N_Xe', 'm_Xe', 'V_gift', 'Sigma_agift', 'eta', 'f', 'k'])  # Write header
    writer.writerows(data)

print(f"Data has been saved to '{file_name}'")


