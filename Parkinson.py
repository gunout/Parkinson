import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Circle, Rectangle
import warnings
warnings.filterwarnings('ignore')

class ParkinsonGenomeAnalyzer:
    def __init__(self):
        self.colors = {'A': '#FF6B6B', 'T': '#4ECDC4', 'C': '#45B7D1', 'G': '#F9A602'}
        self.genes_parkinson = ['SNCA', 'LRRK2', 'PARK2', 'PINK1', 'DJ1', 'GBA', 'VPS35', 'ATP13A2']
        
    def generate_genomic_data(self):
        """Génère des données génomiques simulées pour la maladie de Parkinson"""
        print("🧬 Génération des données génomiques de la maladie de Parkinson...")
        
        # Mutations communes dans la maladie de Parkinson
        mutations_data = {
            'Gene': self.genes_parkinson + ['UCHL1', 'HTRA2', 'PLA2G6', 'FBXO7', 'DNAJC6', 'SYNJ1'],
            'Mutation_Frequency': [0.15, 0.12, 0.08, 0.06, 0.05, 0.10, 0.04, 0.03, 
                                  0.02, 0.02, 0.03, 0.02, 0.02, 0.01],
            'Mutation_Type': ['Missense', 'Missense', 'Deletion', 'Missense', 'Deletion', 
                            'Missense', 'Missense', 'Missense', 'Missense', 'Missense',
                            'Missense', 'Missense', 'Missense', 'Frameshift'],
            'Clinical_Significance': ['Pathogenic', 'Pathogenic', 'Pathogenic', 'Pathogenic',
                                    'Pathogenic', 'Risk_Factor', 'Pathogenic', 'Pathogenic',
                                    'Pathogenic', 'Pathogenic', 'Pathogenic', 'Pathogenic',
                                    'Pathogenic', 'Pathogenic']
        }
        
        return pd.DataFrame(mutations_data)
    
    def create_ctag_diagram(self, df):
        """Crée un diagramme CTAG pour visualiser les mutations"""
        plt.style.use('seaborn-v0_8')
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(20, 16))
        
        # 1. Fréquence des mutations par gène
        self._plot_mutation_frequency(df, ax1)
        
        # 2. Répartition des types de mutations
        self._plot_mutation_types(df, ax2)
        
        # 3. Diagramme CTAG sequence
        self._plot_ctag_sequence(ax3)
        
        # 4. Impact clinique
        self._plot_clinical_impact(df, ax4)
        
        plt.suptitle('Analyse Génomique de la Maladie de Parkinson - Diagramme CTAG', 
                    fontsize=16, fontweight='bold')
        plt.tight_layout()
        plt.savefig('parkinson_ctag_diagram.png', dpi=300, bbox_inches='tight')
        plt.show()
        
        # Générer le rapport d'analyse
        self._generate_genomic_report(df)
    
    def _plot_mutation_frequency(self, df, ax):
        """Plot de la fréquence des mutations par gène"""
        colors = ['#FF6B6B' if gene in ['SNCA', 'LRRK2', 'PARK2'] else '#45B7D1' for gene in df['Gene']]
        
        bars = ax.barh(df['Gene'], df['Mutation_Frequency']*100, 
                      color=colors, alpha=0.8)
        
        # Ajouter les valeurs sur les barres
        for bar in bars:
            width = bar.get_width()
            ax.text(width + 0.5, bar.get_y() + bar.get_height()/2, 
                   f'{width:.1f}%', ha='left', va='center', fontweight='bold', fontsize=9)
        
        ax.set_xlabel('Fréquence des Mutations (%)')
        ax.set_title('Fréquence des Mutations par Gène\n dans la Maladie de Parkinson', 
                    fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3, axis='x')
    
    def _plot_mutation_types(self, df, ax):
        """Plot de la répartition des types de mutations"""
        mutation_counts = df['Mutation_Type'].value_counts()
        
        colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#F9A602', '#6A0572', '#2A9D8F']
        wedges, texts, autotexts = ax.pie(mutation_counts.values, 
                                         labels=mutation_counts.index,
                                         colors=colors[:len(mutation_counts)],
                                         autopct='%1.1f%%',
                                         startangle=90)
        
        for autotext in autotexts:
            autotext.set_color('white')
            autotext.set_fontweight('bold')
        
        ax.set_title('Répartition des Types de Mutations', 
                    fontsize=12, fontweight='bold')
    
    def _plot_ctag_sequence(self, ax):
        """Crée un diagramme CTAG sequence"""
        # Séquence simulée avec mutations spécifiques à Parkinson
        sequence = "ATGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC"
        mutations = [8, 15, 22, 30, 42]  # positions de mutations
        
        ax.set_xlim(0, len(sequence))
        ax.set_ylim(0, 2)
        
        # Dessiner la séquence de base
        for i, base in enumerate(sequence):
            color = self.colors[base]
            rect = Rectangle((i, 0.8), 1, 0.4, facecolor=color, alpha=0.7, edgecolor='black')
            ax.add_patch(rect)
            ax.text(i + 0.5, 1, base, ha='center', va='center', 
                   fontweight='bold', fontsize=8)
        
        # Marquer les mutations
        for mut_pos in mutations:
            circle = Circle((mut_pos + 0.5, 1.6), 0.3, 
                          facecolor='#8B4513', alpha=0.8, edgecolor='#654321')  # Marron pour Parkinson
            ax.add_patch(circle)
            ax.text(mut_pos + 0.5, 1.6, 'M', ha='center', va='center', 
                   fontweight='bold', color='white', fontsize=8)
        
        ax.set_title('Diagramme CTAG - Séquence Génomique avec Mutations Parkinson', 
                    fontsize=12, fontweight='bold')
        ax.set_xlabel('Position dans la séquence')
        ax.set_ylabel('')
        ax.set_yticks([])
        ax.grid(True, alpha=0.3)
    
    def _plot_clinical_impact(self, df, ax):
        """Plot de l'impact clinique des mutations"""
        clinical_counts = df['Clinical_Significance'].value_counts()
        
        colors = ['#E76F51', '#2A9D8F', '#F9A602']
        bars = ax.bar(clinical_counts.index, clinical_counts.values, 
                     color=colors[:len(clinical_counts)], alpha=0.8)
        
        # Ajouter les valeurs sur les barres
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2, height + 0.1,
                   f'{height}', ha='center', va='bottom', fontweight='bold')
        
        ax.set_ylabel('Nombre de Gènes')
        ax.set_title('Impact Clinique des Mutations', 
                    fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3, axis='y')
        
        # Rotation des labels x
        plt.setp(ax.get_xticklabels(), rotation=45, ha='right')
    
    def create_advanced_genomic_analysis(self, df):
        """Crée une analyse génomique avancée"""
        fig = plt.figure(figsize=(18, 12))
        
        # 1. Heatmap des mutations
        ax1 = plt.subplot(2, 3, 1)
        self._plot_mutation_heatmap(df, ax1)
        
        # 2. Distribution des fréquences
        ax2 = plt.subplot(2, 3, 2)
        self._plot_frequency_distribution(df, ax2)
        
        # 3. Réseau d'interaction des gènes
        ax3 = plt.subplot(2, 3, 3)
        self._plot_gene_interaction_network(ax3)
        
        # 4. Profil mutationnel
        ax4 = plt.subplot(2, 3, 4)
        self._plot_mutational_signature(ax4)
        
        # 5. Pathway analysis
        ax5 = plt.subplot(2, 3, 5)
        self._plot_pathway_analysis(df, ax5)
        
        # 6. CTAG detailed
        ax6 = plt.subplot(2, 3, 6)
        self._plot_detailed_ctag(ax6)
        
        plt.suptitle('Analyse Génomique Avancée - Maladie de Parkinson', 
                    fontsize=16, fontweight='bold')
        plt.tight_layout()
        plt.savefig('advanced_parkinson_genomics.png', dpi=300, bbox_inches='tight')
        plt.show()
    
    def _plot_mutation_heatmap(self, df, ax):
        """Crée une heatmap des mutations"""
        # Créer une matrice simulée pour la heatmap
        genes = df['Gene'].tolist()
        mutation_types = df['Mutation_Type'].unique()
        
        # Matrice de fréquence simulée
        heatmap_data = np.random.rand(len(genes), len(mutation_types)) * 0.2
        
        # Renforcer certaines associations spécifiques à Parkinson
        for i, gene in enumerate(genes):
            if gene in ['SNCA', 'LRRK2']:
                heatmap_data[i, 0] = 0.8  # Forte association avec missense
            if gene == 'PARK2':
                heatmap_data[i, 1] = 0.7  # Forte association avec deletion
        
        sns.heatmap(heatmap_data, ax=ax, cmap='YlOrRd', 
                   xticklabels=mutation_types, yticklabels=genes,
                   cbar_kws={'label': 'Fréquence relative'})
        
        ax.set_title('Heatmap des Mutations\npar Gène et Type', fontsize=10, fontweight='bold')
        ax.tick_params(axis='x', rotation=45)
        ax.tick_params(axis='y', rotation=0)
    
    def _plot_frequency_distribution(self, df, ax):
        """Plot de la distribution des fréquences de mutations"""
        ax.hist(df['Mutation_Frequency']*100, bins=8, 
               color='#8B4513', alpha=0.7, edgecolor='black')  # Marron pour Parkinson
        
        ax.axvline(df['Mutation_Frequency'].mean()*100, color='red', 
                  linestyle='--', linewidth=2, label=f'Moyenne: {df["Mutation_Frequency"].mean()*100:.1f}%')
        
        ax.set_xlabel('Fréquence des Mutations (%)')
        ax.set_ylabel('Nombre de Gènes')
        ax.set_title('Distribution des Fréquences de Mutations', 
                    fontsize=10, fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)
    
    def _plot_gene_interaction_network(self, ax):
        """Crée un réseau d'interaction des gènes Parkinson"""
        # Positions simulées pour le réseau
        nodes = {
            'SNCA': (3, 4), 'LRRK2': (5, 3), 'PARK2': (2, 2),
            'PINK1': (4, 1), 'DJ1': (1, 3), 'GBA': (5, 5),
            'VPS35': (3, 5), 'ATP13A2': (1, 1)
        }
        
        # Connexions simulées spécifiques à Parkinson
        edges = [('SNCA', 'LRRK2'), ('SNCA', 'PARK2'), ('PARK2', 'PINK1'),
                ('PINK1', 'DJ1'), ('LRRK2', 'GBA'), ('SNCA', 'VPS35'),
                ('PARK2', 'ATP13A2'), ('GBA', 'VPS35')]
        
        # Dessiner les arêtes
        for edge in edges:
            start = nodes[edge[0]]
            end = nodes[edge[1]]
            ax.plot([start[0], end[0]], [start[1], end[1]], 
                   'gray', alpha=0.6, linewidth=2)
        
        # Dessiner les nœuds
        for gene, pos in nodes.items():
            color = '#8B4513' if gene in ['SNCA', 'LRRK2', 'PARK2'] else '#45B7D1'  # Marron pour gènes principaux
            ax.scatter(pos[0], pos[1], s=300, c=color, alpha=0.8, 
                      edgecolors='black', linewidth=2)
            ax.text(pos[0], pos[1], gene, ha='center', va='center', 
                   fontweight='bold', fontsize=8)
        
        ax.set_xlim(0, 6)
        ax.set_ylim(0, 6)
        ax.set_title('Réseau d\'Interaction des Gènes Parkinson', 
                    fontsize=10, fontweight='bold')
        ax.set_xticks([])
        ax.set_yticks([])
    
    def _plot_mutational_signature(self, ax):
        """Plot des signatures mutationnelles spécifiques à Parkinson"""
        # Signatures mutationnelles simulées pour Parkinson
        contexts = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
        signature1 = [0.15, 0.08, 0.35, 0.12, 0.20, 0.10]  # Signature sporadique
        signature2 = [0.25, 0.12, 0.25, 0.15, 0.15, 0.08]   # Signature familiale
        
        x = np.arange(len(contexts))
        width = 0.35
        
        ax.bar(x - width/2, signature1, width, label='Signature Sporadique', 
               color='#8B4513', alpha=0.8)  # Marron
        ax.bar(x + width/2, signature2, width, label='Signature Familiale', 
               color='#D2691E', alpha=0.8)  # Marron clair
        
        ax.set_xlabel('Type de Mutation')
        ax.set_ylabel('Fréquence relative')
        ax.set_title('Signatures Mutationnelles Parkinson', fontsize=10, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(contexts)
        ax.legend()
        ax.grid(True, alpha=0.3, axis='y')
    
    def _plot_pathway_analysis(self, df, ax):
        """Analyse des voies de signalisation spécifiques à Parkinson"""
        pathways = {
            'Protéostasie': ['SNCA', 'PARK2', 'UCHL1'],
            'Fonction Mitochondriale': ['PINK1', 'DJ1', 'ATP13A2'],
            'Trafic Vésiculaire': ['LRRK2', 'VPS35', 'DNAJC6'],
            'Lysosomal': ['GBA', 'ATP13A2']
        }
        
        pathway_counts = {pathway: len(genes) for pathway, genes in pathways.items()}
        
        colors = ['#8B4513', '#D2691E', '#A0522D', '#CD853F']  # Tons marron
        bars = ax.bar(pathway_counts.keys(), pathway_counts.values(), 
                     color=colors, alpha=0.8)
        
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2, height + 0.1,
                   f'{height}', ha='center', va='bottom', fontweight='bold')
        
        ax.set_ylabel('Nombre de Gènes')
        ax.set_title('Voies de Signalisation Parkinson', 
                    fontsize=10, fontweight='bold')
        ax.tick_params(axis='x', rotation=45)
        ax.grid(True, alpha=0.3, axis='y')
    
    def _plot_detailed_ctag(self, ax):
        """Diagramme CTAG détaillé avec transitions/transversions"""
        # Données simulées pour les transitions/transversions (spécifique Parkinson)
        mutation_types = ['C>T', 'C>G', 'C>A', 'T>C', 'T>G', 'T>A']
        counts = [120, 30, 40, 90, 25, 35]  # Comptes simulés pour Parkinson
        
        colors = ['#8B4513', '#D2691E', '#A0522D', '#CD853F', '#DEB887', '#F4A460']  # Tons marron
        bars = ax.bar(mutation_types, counts, color=colors, alpha=0.8)
        
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2, height + 2,
                   f'{height}', ha='center', va='bottom', fontweight='bold')
        
        ax.set_ylabel('Nombre de Mutations')
        ax.set_title('Spectre des Mutations CTAG Parkinson', 
                    fontsize=10, fontweight='bold')
        ax.grid(True, alpha=0.3, axis='y')
    
    def _generate_genomic_report(self, df):
        """Génère un rapport d'analyse génomique pour Parkinson"""
        print("\n🧬 RAPPORT D'ANALYSE GÉNOMIQUE - MALADIE DE PARKINSON")
        print("=" * 60)
        
        # 1. Gènes les plus fréquemment mutés
        print("\n1. 📊 GÈNES LES PLUS MUTÉS:")
        top_genes = df.nlargest(5, 'Mutation_Frequency')
        for _, row in top_genes.iterrows():
            print(f"   • {row['Gene']}: {row['Mutation_Frequency']*100:.1f}% "
                  f"({row['Mutation_Type']}) - {row['Clinical_Significance']}")
        
        # 2. Analyse des gènes principaux
        print("\n2. 🔍 ANALYSE DES GÈNES PRINCIPAUX:")
        main_genes = df[df['Gene'].isin(['SNCA', 'LRRK2', 'PARK2', 'GBA'])]
        for _, row in main_genes.iterrows():
            print(f"   • {row['Gene']}: Fréquence {row['Mutation_Frequency']*100:.1f}% - {row['Clinical_Significance']}")
        
        # 3. Implications thérapeutiques
        print("\n3. 💊 IMPLICATIONS THÉRAPEUTIQUES:")
        print("   • Thérapies ciblant l'alpha-synucléine: Pour mutations SNCA")
        print("   • Inhibiteurs de LRRK2: En développement clinique")
        print("   • Thérapies enzymatiques: Pour mutations GBA")
        print("   • Antioxydants et protecteurs mitochondriaux")
        
        # 4. Recommandations de dépistage
        print("\n4. 🎯 RECOMMANDATIONS DE DÉPISTAGE:")
        print("   • Séquençage panel Parkinson: SNCA, LRRK2, PARK2, GBA")
        print("   • Conseil génétique: Pour formes familiales")
        print("   • Test GBA: Important pour la réponse thérapeutique")
        
        # 5. Statistiques globales
        print("\n5. 📈 STATISTIQUES GLOBALES:")
        print(f"   • Nombre total de gènes analysés: {len(df)}")
        print(f"   • Fréquence moyenne des mutations: {df['Mutation_Frequency'].mean()*100:.1f}%")
        print(f"   • Gènes pathogènes: {len(df[df['Clinical_Significance'] == 'Pathogenic'])}")
        print(f"   • Facteurs de risque: {len(df[df['Clinical_Significance'] == 'Risk_Factor'])}")

def main():
    """Fonction principale"""
    print("🧬 ANALYSE GÉNOMIQUE DE LA MALADIE DE PARKINSON - DIAGRAMME CTAG")
    print("=" * 60)
    
    # Initialiser l'analyseur
    analyzer = ParkinsonGenomeAnalyzer()
    
    # Générer les données
    genomic_data = analyzer.generate_genomic_data()
    
    # Aperçu des données
    print("\n👀 Aperçu des données génomiques:")
    print(genomic_data.head(8))
    
    # Créer le diagramme CTAG principal
    print("\n📊 Création du diagramme CTAG...")
    analyzer.create_ctag_diagram(genomic_data)
    
    # Créer l'analyse avancée
    print("\n🔬 Création de l'analyse génomique avancée...")
    analyzer.create_advanced_genomic_analysis(genomic_data)
    
    print("\n✅ Analyse génomique terminée!")
    print("📁 Fichiers générés:")
    print("   • parkinson_ctag_diagram.png")
    print("   • advanced_parkinson_genomics.png")

if __name__ == "__main__":
    main()