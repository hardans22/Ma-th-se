import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# Tes données
data = {
    "Inst": list(range(1,25)),
    "FH": [23,8,0,11,0,0,0,0, 12,8,5,56,2,0,1,0, 65,147,209,184,5,11,6,43],
    "FC": [21,20,11,6,16,19,14,22, 32,37,22,28,32,36,31,21, 40,68,24,36,51,45,36,34],
    "DY": [0,0,0,1,0,1,1,2, 2,3,1,3,0,4,1,2, 1,5,1,4,3,1,3,6],
}
df = pd.DataFrame(data)
df.set_index("Inst", inplace=True)

# Définir les blocs de 8
blocks = [(0,8), (8,16), (16,24)]

# Créer une figure avec GridSpec pour contrôler la disposition
fig = plt.figure(figsize=(15, 10))
gs = gridspec.GridSpec(2, 4, figure=fig)  # 2 lignes, 4 colonnes

# Premier graphique (en haut à gauche) - colonnes 0-1
ax1 = fig.add_subplot(gs[0, 0:2])
df_block = df.iloc[blocks[0][0]:blocks[0][1]]
im1 = ax1.imshow(df_block.T, cmap="Blues", aspect="auto")

# Axes
ax1.set_xticks(range(len(df_block.index)))
ax1.set_xticklabels(df_block.index)
ax1.set_yticks(range(len(df_block.columns)))
ax1.set_yticklabels(df_block.columns)

# Ajouter les valeurs
for i, inst in enumerate(df_block.index):
    for j, var in enumerate(df_block.columns):
        ax1.text(i, j, df_block.loc[inst, var], ha="center", va="center", color="black")

ax1.set_title("Test case 1")
ax1.set_xlabel("Instances")
ax1.set_ylabel("Maintenance indicators")
fig.colorbar(im1, ax=ax1, orientation="vertical", fraction=0.05, pad=0.05, label="Values")

# Deuxième graphique (en haut à droite) - colonnes 2-3
ax2 = fig.add_subplot(gs[0, 2:4])
df_block = df.iloc[blocks[1][0]:blocks[1][1]]
im2 = ax2.imshow(df_block.T, cmap="Blues", aspect="auto")

# Axes
ax2.set_xticks(range(len(df_block.index)))
ax2.set_xticklabels(df_block.index)
ax2.set_yticks(range(len(df_block.columns)))
ax2.set_yticklabels(df_block.columns)

# Ajouter les valeurs
for i, inst in enumerate(df_block.index):
    for j, var in enumerate(df_block.columns):
        ax2.text(i, j, df_block.loc[inst, var], ha="center", va="center", color="black")

ax2.set_title("Test case 2")
ax2.set_xlabel("Instances")
ax2.set_ylabel("Maintenance indicators")
fig.colorbar(im2, ax=ax2, orientation="vertical", fraction=0.05, pad=0.05, label="Values")

# Troisième graphique (centré en bas) - colonnes 1-2 (centré)
ax3 = fig.add_subplot(gs[1, 1:3])
df_block = df.iloc[blocks[2][0]:blocks[2][1]]
im3 = ax3.imshow(df_block.T, cmap="Blues", aspect="auto")

# Axes
ax3.set_xticks(range(len(df_block.index)))
ax3.set_xticklabels(df_block.index)
ax3.set_yticks(range(len(df_block.columns)))
ax3.set_yticklabels(df_block.columns)

# Ajouter les valeurs
for i, inst in enumerate(df_block.index):
    for j, var in enumerate(df_block.columns):
        ax3.text(i, j, df_block.loc[inst, var], ha="center", va="center", color="black")

ax3.set_title("Test case 3")
ax3.set_xlabel("Instances")
ax3.set_ylabel("Maintenance indicators")
fig.colorbar(im3, ax=ax3, orientation="vertical", fraction=0.05, pad=0.05, label="Values")

plt.tight_layout()
plt.show()