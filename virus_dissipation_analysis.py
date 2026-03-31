"""
Virus-Specific UV Dissipation Analysis v2

Two detectability scores:
  1. Direct absorption (UVB context): ε280/nt × G-density × genome length
     Relevant to germicidal UV, lab assays with UVB sources.
  2. Photosensitised ROS pathway (UVA/solar context): G-density × genome length
     Under UVA-dominated solar spectra, photon capture occurs via organic
     photosensitisers (humic-like substances, brown carbon) in the microlayer
     matrix, not by nucleic acids directly. All virions in the same matrix
     receive similar ROS flux; the differential is target availability
     (guanine sites for 8-oxoG formation) and total observation sites
     (genome length for phylodynamic detection).

This corrects the spectral alignment issue: ground-level solar UV is
UVA-dominated where nucleic acid absorption is weak.

Adem S. M. Harman, 2026
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ============================================================================
# 1. NUCLEOTIDE ABSORPTION COEFFICIENTS
# ============================================================================
EPS_260 = {'A': 15020, 'U': 9660, 'G': 12080, 'C': 7070}
EPS_280 = {'A': 2020, 'U': 1420, 'G': 8110, 'C': 1520}

# ============================================================================
# 2. VIRUS GENOMES
# ============================================================================
VIRUSES = [
    ("SARS-CoV-2",       29903, 29.9, 32.1, 19.6, 18.4, 1, "Coronaviridae"),
    ("SARS-CoV-1",       29751, 28.8, 32.5, 20.7, 18.0, 1, "Coronaviridae"),
    ("MERS-CoV",         30119, 27.4, 31.8, 21.0, 19.8, 1, "Coronaviridae"),
    ("HCoV-229E",        27317, 28.7, 32.5, 20.5, 18.3, 1, "Coronaviridae"),
    ("HCoV-OC43",        30741, 28.7, 33.2, 20.2, 17.9, 1, "Coronaviridae"),
    ("Influenza A",      13588, 33.1, 26.9, 22.8, 17.2, 8, "Orthomyxoviridae"),
    ("Influenza B",      14452, 33.0, 28.0, 22.0, 17.0, 8, "Orthomyxoviridae"),
    ("RSV",              15222, 32.3, 27.5, 20.3, 19.9, 1, "Pneumoviridae"),
    ("Measles",          15894, 33.5, 26.2, 23.3, 17.0, 1, "Paramyxoviridae"),
    ("Mumps",            15384, 30.0, 24.0, 26.0, 20.0, 1, "Paramyxoviridae"),
    ("Ebola (Zaire)",    18959, 32.0, 27.0, 21.0, 20.0, 1, "Filoviridae"),
    ("Dengue-2",         10723, 30.0, 23.0, 26.0, 21.0, 1, "Flaviviridae"),
    ("Zika",             10794, 28.5, 22.5, 27.0, 22.0, 1, "Flaviviridae"),
    ("Hepatitis C",       9646, 27.0, 21.0, 29.0, 23.0, 1, "Flaviviridae"),
    ("Rhinovirus A",      7152, 28.5, 26.0, 24.0, 21.5, 1, "Picornaviridae"),
    ("Enterovirus D68",   7367, 28.0, 24.0, 26.0, 22.0, 1, "Picornaviridae"),
    ("Norovirus GII",     7654, 27.0, 25.0, 26.0, 22.0, 1, "Caliciviridae"),
    ("Rabies",           11932, 28.5, 24.5, 26.0, 21.0, 1, "Rhabdoviridae"),
    ("HIV-1",             9719, 35.8, 22.4, 24.0, 17.8, 1, "Retroviridae"),
]

# ============================================================================
# 3. COMPUTE BOTH SCORES
# ============================================================================
results = []

for name, length, pA, pU, pG, pC, segments, family in VIRUSES:
    fA, fU, fG, fC = pA/100, pU/100, pG/100, pC/100
    
    eps_280 = fA*EPS_280['A'] + fU*EPS_280['U'] + fG*EPS_280['G'] + fC*EPS_280['C']
    
    g_content = fG
    g_per_kb = g_content * 1000
    
    # Score 1: Direct absorption (UVB/lab context)
    # Photon capture by RNA × target sites × observation sites
    score_direct = eps_280 * g_content * length / 1e9
    
    # Score 2: Photosensitised ROS (UVA/solar context)
    # Matrix delivers ROS uniformly; differential = target × observation
    # G sites × genome length (no ε term — photon capture is matrix-driven)
    score_ros = g_content * length / 1e3
    
    # Score 3: Per-genome G-site target density (lab sensitivity)
    # Total G sites = target count for 8-oxoG
    total_g = g_content * length
    
    results.append({
        'name': name, 'family': family, 'length': length,
        'segments': segments, 'g_content': g_content, 'g_per_kb': g_per_kb,
        'eps_280': eps_280, 'total_g': total_g,
        'score_direct': score_direct,
        'score_ros': score_ros,
    })

# ============================================================================
# 4. RESULTS
# ============================================================================
print("=" * 85)
print("VIRUS-SPECIFIC UV DISSIPATION ANALYSIS v2")
print("Dual scoring: Direct UVB absorption vs Photosensitised UVA/ROS pathway")
print("=" * 85)

# Sort by ROS score (the spectrally correct one for solar exposure)
results.sort(key=lambda x: x['score_ros'], reverse=True)

print(f"\n{'Virus':<18} {'Length':>7} {'G%':>5} {'G/kb':>5} {'ε₂₈₀/nt':>8} "
      f"{'Total G':>7} {'Direct':>8} {'ROS':>8}")
print("-" * 85)
for r in results:
    print(f"{r['name']:<18} {r['length']:>7,} {r['g_content']*100:>4.1f}% "
          f"{r['g_per_kb']:>5.0f} {r['eps_280']:>8,.0f} "
          f"{r['total_g']:>7,.0f} {r['score_direct']:>8.1f} {r['score_ros']:>8.2f}")

# ============================================================================
# 5. RANKING COMPARISON
# ============================================================================
print("\n" + "=" * 85)
print("RANKING COMPARISON: Does spectral correction change priorities?")
print("=" * 85)

results_direct = sorted(results, key=lambda x: x['score_direct'], reverse=True)
results_ros = sorted(results, key=lambda x: x['score_ros'], reverse=True)

print(f"\n{'Rank':>4}  {'Direct (UVB)':^22}  {'Photosensitised (UVA/ROS)':^26}")
print("-" * 60)
for i in range(len(results)):
    d = results_direct[i]
    r = results_ros[i]
    match = "✓" if d['name'] == r['name'] else ""
    print(f"{i+1:>4}  {d['name']:<22}  {r['name']:<26} {match}")

# Count rank changes
direct_order = [r['name'] for r in results_direct]
ros_order = [r['name'] for r in results_ros]

rank_shifts = []
for name in direct_order:
    d_rank = direct_order.index(name)
    r_rank = ros_order.index(name)
    shift = abs(d_rank - r_rank)
    if shift > 0:
        rank_shifts.append((name, d_rank+1, r_rank+1, shift))

print(f"\nRank changes > 0:")
for name, dr, rr, shift in sorted(rank_shifts, key=lambda x: -x[3]):
    direction = "↑" if rr < dr else "↓"
    print(f"  {name:<20} Direct: #{dr} → ROS: #{rr}  ({direction}{shift})")

# Key comparison
corona_direct = np.mean([r['score_direct'] for r in results if r['family'] == 'Coronaviridae'])
corona_ros = np.mean([r['score_ros'] for r in results if r['family'] == 'Coronaviridae'])
other_direct = np.mean([r['score_direct'] for r in results if r['family'] != 'Coronaviridae'])
other_ros = np.mean([r['score_ros'] for r in results if r['family'] != 'Coronaviridae'])

print(f"\nCoronavirus advantage:")
print(f"  Direct (UVB):  {corona_direct/other_direct:.1f}× non-coronavirus mean")
print(f"  ROS (UVA):     {corona_ros/other_ros:.1f}× non-coronavirus mean")

flavi_direct = np.mean([r['score_direct'] for r in results if r['family'] == 'Flaviviridae'])
flavi_ros = np.mean([r['score_ros'] for r in results if r['family'] == 'Flaviviridae'])
print(f"\nFlaviviridae:")
print(f"  Direct rank advantage from high ε₂₈₀: YES (high per-nt absorption)")
print(f"  ROS rank: per-nt absorption irrelevant; drops to target density only")
print(f"  Flavivirus mean score - Direct: {flavi_direct:.1f}, ROS: {flavi_ros:.2f}")

# ============================================================================
# 6. KEY FINDING
# ============================================================================
print(f"\n" + "=" * 85)
print("KEY FINDING")
print("=" * 85)
print("""
The spectral correction PRESERVES the top-level finding:

  Coronaviruses remain the optimal phylodynamic targets under BOTH
  scoring models, because genome length and total G-site count dominate
  regardless of whether photon capture occurs via direct RNA absorption
  (UVB) or photosensitised ROS generation (UVA).

What changes:

  1. The per-nucleotide ε₂₈₀ axis (Panel B in the original figure) becomes
     LESS relevant under solar/UVA conditions — it matters for laboratory
     UVB assays but not for environmental UVA-driven oxidative damage.

  2. Flaviviruses lose their per-nucleotide advantage under the ROS model
     because all virions in the same matrix receive similar ROS flux. Their
     high G content still makes them good lab targets but the advantage is
     smaller and driven purely by target density, not photon capture.

  3. The DIFFERENTIAL PREDICTION sharpens:
     - Lab (UVB): both photon capture AND target density vary across viruses
       → flaviviruses show strongest per-genome effects
     - Field (UVA/solar): photon capture is matrix-driven (constant across
       viruses in same environment) → only target density and genome length
       matter → coronaviruses dominate even more strongly

  4. The microlayer argument is STRENGTHENED: under solar conditions, the
     organic composition of the microlayer (photosensitiser content) drives
     the photon-to-ROS conversion, making the microlayer chemistry the
     central variable — exactly as dissipation theory predicts.
""")

# ============================================================================
# 7. FIGURE
# ============================================================================
family_colors = {
    'Coronaviridae': '#c0392b', 'Orthomyxoviridae': '#2980b9',
    'Paramyxoviridae': '#27ae60', 'Pneumoviridae': '#27ae60',
    'Filoviridae': '#8e44ad', 'Flaviviridae': '#e67e22',
    'Picornaviridae': '#7f8c8d', 'Caliciviridae': '#95a5a6',
    'Rhabdoviridae': '#34495e', 'Retroviridae': '#d35400',
}

fig, axes = plt.subplots(1, 3, figsize=(15, 6))

# Panel A: Genome length vs total G sites
ax = axes[0]
for r in results:
    c = family_colors.get(r['family'], '#95a5a6')
    ax.scatter(r['length']/1000, r['total_g'], c=c, s=80, 
              edgecolors='white', linewidth=0.5, zorder=3)
    if r['name'] in ['SARS-CoV-2', 'Influenza A', 'Ebola (Zaire)',
                      'Hepatitis C', 'Rhinovirus A', 'HIV-1', 'Zika']:
        ax.annotate(r['name'], (r['length']/1000, r['total_g']),
                   fontsize=7, ha='left', va='bottom',
                   xytext=(4, 4), textcoords='offset points')
ax.set_xlabel('Genome length (kb)', fontsize=10)
ax.set_ylabel('Total guanine sites\n(8-oxoG targets)', fontsize=10)
ax.set_title('A.  Target availability', fontsize=11, fontweight='bold', loc='left')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel B: Direct vs ROS detectability (rank comparison)
ax = axes[1]
for r in results:
    c = family_colors.get(r['family'], '#95a5a6')
    ax.scatter(r['score_direct'], r['score_ros'], c=c, s=80,
              edgecolors='white', linewidth=0.5, zorder=3)
    if r['name'] in ['SARS-CoV-2', 'MERS-CoV', 'Hepatitis C', 
                      'Influenza A', 'Zika', 'Rhinovirus A', 'HIV-1']:
        ax.annotate(r['name'], (r['score_direct'], r['score_ros']),
                   fontsize=7, ha='left', va='bottom',
                   xytext=(4, 4), textcoords='offset points')

# Add diagonal reference
max_d = max(r['score_direct'] for r in results)
max_r = max(r['score_ros'] for r in results)
ax.set_xlabel('Direct UVB score', fontsize=10)
ax.set_ylabel('Photosensitised UVA/ROS score', fontsize=10)
ax.set_title('B.  UVB vs UVA pathway scores', fontsize=11, fontweight='bold', loc='left')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel C: Dual ranking
ax = axes[2]
names_ros = [r['name'] for r in results_ros]
scores_ros = [r['score_ros'] for r in results_ros]
colors_ros = [family_colors.get(r['family'], '#95a5a6') for r in results_ros]

bars = ax.barh(range(len(names_ros)), scores_ros, color=colors_ros, 
               edgecolor='white', linewidth=0.5, alpha=0.8)
ax.set_yticks(range(len(names_ros)))
ax.set_yticklabels(names_ros, fontsize=8)
ax.set_xlabel('Photosensitised ROS score\n(solar/environmental context)', fontsize=10)
ax.set_title('C.  H3 detectability (solar UV)', fontsize=11, fontweight='bold', loc='left')
ax.invert_yaxis()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

from matplotlib.patches import Patch
families_shown = ['Coronaviridae', 'Orthomyxoviridae', 'Paramyxoviridae',
                  'Filoviridae', 'Flaviviridae', 'Picornaviridae', 'Retroviridae']
legend_elements = [Patch(facecolor=family_colors[f], label=f) for f in families_shown]
fig.legend(handles=legend_elements, loc='lower center', ncol=4, fontsize=8,
          bbox_to_anchor=(0.5, -0.02))

plt.tight_layout()
plt.subplots_adjust(bottom=0.12)
plt.savefig('/home/claude/virus_dissipation_figure_v2.pdf', dpi=300, bbox_inches='tight')
plt.savefig('/home/claude/virus_dissipation_figure_v2.png', dpi=200, bbox_inches='tight')
print("Figure saved.")
