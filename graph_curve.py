import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# === Données ===
data = {
"Instance": list(range(1,17)),
# Minimize legal remaining FH -> FC, DY
"FH_group_FC": [18,10,11,6,16,19,14,22,19,16,14,17,28,31,23,21],
"FH_group_DY": [0,0,0,1,0,1,1,2,2,3,1,1,0,4,1,1],
# Minimize legal remaining FC -> FH, DY
"FC_group_FH": [794,286,311,382,328,321,218,534,255,367,338,392,946,741,247,930],
"FC_group_DY": [0,0,0,1,0,0,0,1,2,3,0,0,0,0,0,0],
# Minimize legal remaining DY -> FH, FC
"DY_group_FH": [445,279,261,84,2278,716,195,231,436,879,601,1056,2009,1063,711,725],
"DY_group_FC": [17,9,14,7,31,14,10,7,14,10,15,13,45,27,18,22],
}
df = pd.DataFrame(data)

# === 1. Minimize FH ===
fig1, ax1 = plt.subplots(figsize=(8,4.5))
ax1.plot(df["Instance"], df["FH_group_FC"], color="tab:blue", marker='o', label="FC (remaining)")
ax1.set_xlabel("Instance")
ax1.set_ylabel("FC (remaining)")
ax1.set_title("Minimize legal remaining FH")
ax1.grid(True, linestyle=':', linewidth=0.5)
ax1.set_xlim(0.8, 16.2)
ax1.set_xticks(range(1, 17))

# Définir les limites de l'axe gauche avec un petit espace
ax1_min = -max(df["FH_group_FC"]) * 0.02  # 2% d'espace en bas
ax1_max = max(df["FH_group_FC"]) * 1.0   # 5% d'espace en haut
ax1.set_ylim(ax1_min, ax1_max)

# Créer l'axe de droite avec cohérence (ratio environ 10:1)
ax1b = ax1.twinx()
ax1b.plot(df["Instance"], df["FH_group_DY"], color="tab:orange", marker='o', label="DY (remaining)")
ax1b.set_ylabel("DY (remaining)")

# Calculer les limites de l'axe droit pour avoir une cohérence avec espacement
ratio = ax1_max / max(max(df["FH_group_DY"]), 1) * 0.1  # Ajuster le ratio
ax1b_max = ax1_max / 10  # 10:1 ratio approximatif
ax1b_min = -ax1b_max * 0.02  # Petit espace en bas proportionnel
ax1b.set_ylim(ax1b_min, ax1b_max * 1.05)  # Petit espace en haut

fig1.tight_layout()

# === 2. Minimize FC ===
fig2, ax2 = plt.subplots(figsize=(8,4.5))
ax2.plot(df["Instance"], df["FC_group_FH"], color="tab:blue", marker='o', label="FH (remaining)")
ax2.set_xlabel("Instance")
ax2.set_ylabel("FH (remaining)")
ax2.set_title("Minimize legal remaining FC")
ax2.grid(True, linestyle=':', linewidth=0.5)
ax2.set_xlim(0.8, 16.2)
ax2.set_xticks(range(1, 17))

# Définir les limites de l'axe gauche avec un petit espace
ax2_min = -max(df["FC_group_FH"]) * 0.02  # 2% d'espace en bas
ax2_max = max(df["FC_group_FH"]) * 1.05   # 5% d'espace en haut
ax2.set_ylim(ax2_min, ax2_max)

# Créer l'axe de droite avec cohérence
ax2b = ax2.twinx()
ax2b.plot(df["Instance"], df["FC_group_DY"], color="tab:orange", marker='o', label="DY (remaining)")
ax2b.set_ylabel("DY (remaining)")

# Calculer les limites pour avoir une cohérence visuelle avec espacement
ax2b_max = ax2_max / 300  # Ratio approximatif basé sur les données
ax2b_min = -ax2b_max * 0.02  # Petit espace en bas proportionnel
ax2b.set_ylim(ax2b_min, ax2b_max * 1.05)  # Petit espace en haut

fig2.tight_layout()

# === 3. Minimize DY ===
fig3, ax3 = plt.subplots(figsize=(8,4.5))
ax3.plot(df["Instance"], df["DY_group_FH"], color="tab:blue", marker='o', label="FH (remaining)")
ax3.set_xlabel("Instance")
ax3.set_ylabel("FH (remaining)")
ax3.set_title("Minimize legal remaining DY")
ax3.grid(True, linestyle=':', linewidth=0.5)
ax3.set_xlim(0.8, 16.2)
ax3.set_xticks(range(1, 17))

# Définir les limites de l'axe gauche avec un petit espace
ax3_min = -max(df["DY_group_FH"]) * 0.02  # 2% d'espace en bas
ax3_max = max(df["DY_group_FH"]) * 1.05   # 5% d'espace en haut
ax3.set_ylim(ax3_min, ax3_max)

# Créer l'axe de droite avec cohérence
ax3b = ax3.twinx()
ax3b.plot(df["Instance"], df["DY_group_FC"], color="tab:orange", marker='o', label="FC (remaining)")
ax3b.set_ylabel("FC (remaining)")

# Calculer les limites pour avoir une cohérence visuelle avec espacement
ax3b_max = ax3_max / 50  # Ratio approximatif basé sur les données
ax3b_min = -ax3b_max * 0.02  # Petit espace en bas proportionnel
ax3b.set_ylim(ax3b_min, ax3b_max * 1.05)  # Petit espace en haut

fig3.tight_layout()

plt.show()

# Afficher les ratios utilisés pour information
print(f"Graphique 1 - Ratio FH_group_FC/FH_group_DY: {ax1_max/ax1b_max:.1f}")
print(f"Graphique 2 - Ratio FC_group_FH/FC_group_DY: {ax2_max/ax2b_max:.1f}")
print(f"Graphique 3 - Ratio DY_group_FH/DY_group_FC: {ax3_max/ax3b_max:.1f}")