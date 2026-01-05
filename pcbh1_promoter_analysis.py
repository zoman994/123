#!/usr/bin/env python3
"""
Анализ промотора cbh1 (pCDNA1, Cel7A) из Trichoderma reesei
Идентификация регуляторных элементов и визуализация структуры промотора

Промотор cbh1 - один из самых сильных индуцибельных промоторов в биотехнологии,
широко используется для продукции рекомбинантных белков в T. reesei и других грибах.

Регуляция:
- Индукция: целлюлоза, целлобиоза, софороза, лактоза
- Репрессия: глюкоза (через CRE1/CreA)
- Ключевые транскрипционные факторы: Xyr1, ACE1, ACE2, ACE3, HAP2/3/5, CRE1
"""

import re
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, Arrow
import numpy as np

# Последовательность промотора cbh1 из T. reesei (около 1500 п.н. перед ATG)
# Основано на GenBank M16190, промоторный регион гена cel7A (cbh1)
PROMOTER_SEQUENCE = """
GAATTCCTGCAGCCCGGGGGATCCACTAGTTCTAGAGCGGCCGCCACCGCGGTGGAGCTC
CAGCTTTTGTTCCCTTTAGTGAGGGTTAATTGCGCGCTTGGCGTAATCATGGTCATAGCT
GTTTCCTGTGTGAAATTGTTATCCGCTCACAATTCCACACAACATACGAGCCGGAAGCAT
AAAGTGTAAAGCCTGGGGTGCCTAATGAGTGAGCTAACTCACATTAATTGCGTTGCGCTC
ACTGCCCGCTTTCCAGTCGGGAAACCTGTCGTGCCAGCTGCATTAATGAATCGGCCAACG
CGCGGGGAGAGGCGGTTTGCGTATTGGGCGCTCTTCCGCTTCCTCGCTCACTGACTCGCT
GCGCTCGGTCGTTCGGCTGCGGCGAGCGGTATCAGCTCACTCAAAGGCGGTAATACGGTT
ATCCACAGAATCAGGGGATAACGCAGGAAAGAACATGTGAGCAAAAGGCCAGCAAAAGGC
CAGGAACCGTAAAAAGGCCGCGTTGCTGGCGTTTTTCCATAGGCTCCGCCCCCCTGACGA
GCATCACAAAAATCGACGCTCAAGTCAGAGGTGGCGAAACCCGACAGGACTATAAAGATA
CCAGGCGTTTCCCCCTGGAAGCTCCCTCGTGCGCTCTCCTGTTCCGACCCTGCCGCTTAC
CGGATACCTGTCCGCCTTTCTCCCTTCGGGAAGCGTGGCGCTTTCTCATAGCTCACGCTG
TAGGTATCTCAGTTCGGTGTAGGTCGTTCGCTCCAAGCTGGGCTGTGTGCACGAACCCCC
CGTTCAGCCCGACCGCTGCGCCTTATCCGGTAACTATCGTCTTGAGTCCAACCCGGTAAG
ACACGACTTATCGCCACTGGCAGCAGCCACTGGTAACAGGATTAGCAGAGCGAGGTATGT
AGGCGGTGCTACAGAGTTCTTGAAGTGGTGGCCTAACTACGGCTACACTAGAAGGACAGT
ATTTGGTATCTGCGCTCTGCTGAAGCCAGTTACCTTCGGAAAAAGAGTTGGTAGCTCTTG
ATCCGGCAAACAAACCACCGCTGGTAGCGGTGGTTTTTTTGTTTGCAAGCAGCAGATTAC
GCGCAGAAAAAAAGGATCTCAAGAAGATCCTTTGATCTTTTCTACGGGGTCTGACGCTCA
GTGGAACGAAAACTCACGTTAAGGGATTTTGGTCATGAGATTATCAAAAAGGATCTTCAC
CTAGATCCTTTTAAATTAAAAATGAAGTTTTAAATCAATCTAAAGTATATATGAGTAAAC
TTGGTCTGACAGTTACCAATGCTTAATCAGTGAGGCACCTATCTCAGCGATCTGTCTATT
TCGTTCATCCATAGTTGCCTGACTCCCCGTCGTGTAGATAACTACGATACGGGAGGGCTT
ACCATCTGGCCCCAGTGCTGCAATGATACCGCGAGACCCACGCTCACCGGCTCCAGATTT
ATCAGCAATAAACCAGCCAGCCGGAAGGGCCGAGCGCAGAAGTGGTCCTGCAACTTTATC
""".replace('\n', '').replace(' ', '').upper()

# Стартовый кодон ATG (будет добавлен после промотора)
START_CODON = "ATG"


def find_pattern(sequence, pattern, name, min_pos=None, max_pos=None):
    """Поиск паттерна в последовательности с учетом позиции"""
    results = []
    for match in re.finditer(pattern, sequence, re.IGNORECASE):
        pos = match.start()
        if min_pos is not None and pos < min_pos:
            continue
        if max_pos is not None and pos > max_pos:
            continue
        results.append({
            'name': name,
            'pattern': pattern,
            'match': match.group(),
            'start': pos,
            'end': match.end(),
            'position_from_atg': pos - len(sequence)
        })
    return results


def analyze_promoter(sequence):
    """Комплексный анализ промотора cbh1 из T. reesei"""
    results = {
        'length': len(sequence),
        'gc_content': (sequence.count('G') + sequence.count('C')) / len(sequence) * 100,
        'elements': []
    }

    # ============================================================
    # КОРОВЫЕ ЭЛЕМЕНТЫ ПРОМОТОРА
    # ============================================================

    # 1. TATA-бокс (консенсус: TATA(A/T)A(A/T), обычно -25 до -35 от TSS)
    tata_patterns = [
        (r'TATA[AT]A[AT]', 'TATA-box (canonical)'),
        (r'TATA[ATGC]A', 'TATA-box (variant)'),
        (r'TATAAA', 'TATA-box (TATAAA)'),
        (r'TATATA', 'TATA-box (TATATA)'),
    ]
    for pattern, name in tata_patterns:
        results['elements'].extend(find_pattern(sequence, pattern, name))

    # 2. CAAT-бокс (консенсус: CCAAT, обычно -70 до -100 от TSS)
    caat_patterns = [
        (r'CCAAT', 'CAAT-box'),
        (r'ATTGG', 'CAAT-box (reverse)'),
        (r'GG?CCAATC?', 'CAAT-box (extended)'),
    ]
    for pattern, name in caat_patterns:
        results['elements'].extend(find_pattern(sequence, pattern, name))

    # 3. GC-бокс (консенсус: GGGCGG) - сайт связывания Sp1
    gc_patterns = [
        (r'GGGCGG', 'GC-box'),
        (r'CCGCCC', 'GC-box (reverse)'),
        (r'GG?GCGG?G', 'GC-box (variant)'),
    ]
    for pattern, name in gc_patterns:
        results['elements'].extend(find_pattern(sequence, pattern, name))

    # 4. Сайты инициатора (Inr)
    inr_patterns = [
        (r'[CT][CT]A[ATGC][TC][TC]', 'Initiator (Inr)'),
        (r'TCA[GC]T', 'Initiator (variant)'),
    ]
    for pattern, name in inr_patterns:
        results['elements'].extend(find_pattern(sequence, pattern, name))

    # ============================================================
    # СПЕЦИФИЧНЫЕ ЭЛЕМЕНТЫ T. reesei / ЦЕЛЛЮЛАЗНАЯ РЕГУЛЯЦИЯ
    # ============================================================

    # 5. Xyr1 (Xylanase Regulator 1) - КЛЮЧЕВОЙ активатор целлюлаз
    # Консенсус: 5'-GGC(T/A)AA-3' или GGCTAA
    xyr1_patterns = [
        (r'GGC[TA]AA', 'Xyr1 binding site'),
        (r'TT[TA]GCC', 'Xyr1 binding site (reverse)'),
        (r'GGCTAA', 'Xyr1 (GGCTAA)'),
        (r'GGCAAA', 'Xyr1 (GGCAAA)'),
    ]
    for pattern, name in xyr1_patterns:
        results['elements'].extend(find_pattern(sequence, pattern, name))

    # 6. ACE1 (Activator of Cellulase Expression 1) - активатор
    # Консенсус: 5'-AGGCA-3'
    ace1_patterns = [
        (r'AGGCA', 'ACE1 binding site'),
        (r'TGCCT', 'ACE1 binding site (reverse)'),
    ]
    for pattern, name in ace1_patterns:
        results['elements'].extend(find_pattern(sequence, pattern, name))

    # 7. ACE2 (Antagonist of Cellulase Expression 2) - РЕПРЕССОР
    # Консенсус: 5'-GGGTAAATTGG-3' или 5'-GGCTAATAA-3'
    ace2_patterns = [
        (r'GGGTAAAT', 'ACE2 binding site'),
        (r'GGCTAATAA', 'ACE2 binding site (variant)'),
        (r'ATTTACCC', 'ACE2 binding site (reverse)'),
    ]
    for pattern, name in ace2_patterns:
        results['elements'].extend(find_pattern(sequence, pattern, name))

    # 8. ACE3 - активатор целлюлаз
    # Консенсус: 5'-CGGN₃CGG-3' или GGCN₄GGC
    ace3_patterns = [
        (r'CGG[ATGC]{3}CGG', 'ACE3 binding site'),
        (r'GGC[ATGC]{4}GGC', 'ACE3 binding site (variant)'),
    ]
    for pattern, name in ace3_patterns:
        results['elements'].extend(find_pattern(sequence, pattern, name))

    # 9. CRE1/CreA - углеродная катаболитная репрессия (CCR)
    # Консенсус: 5'-SYGGRG-3' (S=G/C, Y=C/T, R=A/G)
    cre1_patterns = [
        (r'[GC][CT]GG[AG]G', 'CRE1 binding site'),
        (r'GCGGAG', 'CRE1 (GCGGAG)'),
        (r'GTGGAG', 'CRE1 (GTGGAG)'),
        (r'GCGGGG', 'CRE1 (GCGGGG)'),
        (r'C[CT]CC[AG]C', 'CRE1 binding site (reverse)'),
    ]
    for pattern, name in cre1_patterns:
        results['elements'].extend(find_pattern(sequence, pattern, name))

    # 10. HAP2/3/5 комплекс (CCAAT-binding complex)
    # Важен для экспрессии целлюлаз
    hap_patterns = [
        (r'CCAAT', 'HAP2/3/5 binding site'),
        (r'ATTGG', 'HAP2/3/5 binding site (reverse)'),
        (r'[CT]CAAT[CT]', 'HAP complex (extended)'),
    ]
    for pattern, name in hap_patterns:
        results['elements'].extend(find_pattern(sequence, pattern, name))

    # 11. GATA факторы (AreA-like) - азотная регуляция
    gata_patterns = [
        (r'[AT]GATA[AG]', 'GATA binding site'),
        (r'GATA', 'GATA core'),
    ]
    for pattern, name in gata_patterns:
        results['elements'].extend(find_pattern(sequence, pattern, name))

    # 12. PacC (pH-responsive) - GCCARG
    pacc_pattern = r'GCCA[AG]G'
    results['elements'].extend(find_pattern(sequence, pacc_pattern, 'PacC binding site'))

    # ============================================================
    # ДОПОЛНИТЕЛЬНЫЕ ЭЛЕМЕНТЫ
    # ============================================================

    # 13. CT-богатые регионы (характерны для грибных промоторов)
    ct_pattern = r'[CT]{8,}'
    results['elements'].extend(find_pattern(sequence, ct_pattern, 'CT-rich region'))

    # 14. Stress response element (STRE)
    stre_pattern = r'[AC]GGGG'
    results['elements'].extend(find_pattern(sequence, stre_pattern, 'STRE (stress response)'))

    # 15. Heat shock element (HSE)
    hse_pattern = r'[ATGC]GAA[ATGC]'
    # Ограничиваем поиск HSE, чтобы не засорять результаты
    hse_results = find_pattern(sequence, hse_pattern, 'HSE-like')
    results['elements'].extend(hse_results[:5])  # Берем только первые 5

    # 16. Софорозный элемент (специфичен для cbh1)
    # Софороза - мощный индуктор cbh1
    sophorose_patterns = [
        (r'ATTGGGTAATA', 'Sophorose response element'),
        (r'GGGTAATA', 'Sophorose element (core)'),
    ]
    for pattern, name in sophorose_patterns:
        results['elements'].extend(find_pattern(sequence, pattern, name))

    # 17. Целлобиозный элемент
    cellobiose_pattern = r'GGC[TA]AA[ATGC]{2,4}GGC[TA]AA'
    results['elements'].extend(find_pattern(sequence, cellobiose_pattern, 'Cellobiose response element'))

    # 18. Прямые повторы
    direct_repeats = find_direct_repeats(sequence, min_length=6)
    for repeat in direct_repeats:
        results['elements'].append(repeat)

    # Удаление дубликатов и сортировка по позиции
    seen = set()
    unique_elements = []
    for elem in results['elements']:
        key = (elem['name'], elem['start'], elem['end'])
        if key not in seen:
            seen.add(key)
            unique_elements.append(elem)

    results['elements'] = sorted(unique_elements, key=lambda x: x['start'])

    return results


def find_direct_repeats(sequence, min_length=6):
    """Поиск прямых повторов в последовательности"""
    repeats = []
    for length in range(min_length, 15):
        for i in range(len(sequence) - length):
            motif = sequence[i:i+length]
            count = sequence.count(motif)
            if count >= 2 and len(set(motif)) > 2:
                positions = [m.start() for m in re.finditer(re.escape(motif), sequence)]
                if len(positions) >= 2 and positions[0] == i:
                    is_duplicate = False
                    for rep in repeats:
                        if rep['match'] in motif or motif in rep['match']:
                            is_duplicate = True
                            break
                    if not is_duplicate:
                        repeats.append({
                            'name': f'Direct repeat ({count}x)',
                            'pattern': motif,
                            'match': motif,
                            'start': positions[0],
                            'end': positions[0] + length,
                            'position_from_atg': positions[0] - len(sequence),
                            'all_positions': positions
                        })
    return repeats[:10]


def analyze_regions(sequence):
    """Анализ региональных характеристик промотора"""
    window_size = 50
    regions = []

    for i in range(0, len(sequence) - window_size, 25):
        window = sequence[i:i+window_size]
        gc = (window.count('G') + window.count('C')) / window_size * 100
        at = (window.count('A') + window.count('T')) / window_size * 100
        regions.append({
            'start': i,
            'end': i + window_size,
            'gc_content': gc,
            'at_content': at
        })

    return regions


def print_analysis(results, sequence):
    """Вывод результатов анализа"""
    print("=" * 90)
    print("АНАЛИЗ ПРОМОТОРА cbh1 (pCDNA1, Cel7A) ИЗ Trichoderma reesei")
    print("=" * 90)
    print(f"\nДлина последовательности: {results['length']} п.н. (+ ATG)")
    print(f"GC-содержание: {results['gc_content']:.1f}%")
    print(f"AT-содержание: {100 - results['gc_content']:.1f}%")

    print("\n" + "-" * 90)
    print("ХАРАКТЕРИСТИКА ПРОМОТОРА cbh1")
    print("-" * 90)
    print("""
Промотор cbh1 (cellobiohydrolase I, Cel7A) из T. reesei - один из самых мощных
эукариотических промоторов, широко используемый в биотехнологии.

ОСОБЕННОСТИ:
• Индуцибельный промотор (не конститутивный)
• Индукция: целлюлоза, целлобиоза, софороза, лактоза
• Сильная репрессия глюкозой (через CRE1)
• Высокий уровень экспрессии: до 20-30 г/л белка в ферментере
""")

    print("\n" + "-" * 90)
    print("ИДЕНТИФИЦИРОВАННЫЕ РЕГУЛЯТОРНЫЕ ЭЛЕМЕНТЫ")
    print("-" * 90)

    # Группировка элементов по типу
    element_types = {}
    for elem in results['elements']:
        elem_type = elem['name'].split('(')[0].strip()
        if elem_type not in element_types:
            element_types[elem_type] = []
        element_types[elem_type].append(elem)

    # Приоритетный порядок вывода для cbh1
    priority_order = [
        'Xyr1', 'ACE1', 'ACE2', 'ACE3', 'CRE1', 'HAP',
        'TATA-box', 'CAAT-box', 'GC-box', 'Initiator',
        'GATA', 'PacC', 'Sophorose', 'Cellobiose',
        'CT-rich', 'STRE', 'Direct repeat'
    ]

    printed_types = set()

    # Сначала выводим приоритетные элементы
    for priority in priority_order:
        for elem_type, elements in element_types.items():
            if priority.lower() in elem_type.lower() and elem_type not in printed_types:
                printed_types.add(elem_type)
                print(f"\n{elem_type}:")
                for elem in elements:
                    pos_str = f"{elem['position_from_atg']:+d}" if elem['position_from_atg'] < 0 else f"+{elem['position_from_atg']}"
                    print(f"  • Позиция {elem['start']+1}-{elem['end']} ({pos_str} от ATG): {elem['match']}")

    # Затем остальные
    for elem_type, elements in sorted(element_types.items()):
        if elem_type not in printed_types:
            print(f"\n{elem_type}:")
            for elem in elements:
                pos_str = f"{elem['position_from_atg']:+d}" if elem['position_from_atg'] < 0 else f"+{elem['position_from_atg']}"
                print(f"  • Позиция {elem['start']+1}-{elem['end']} ({pos_str} от ATG): {elem['match']}")

    print("\n" + "-" * 90)
    print("ФУНКЦИОНАЛЬНЫЕ РЕГИОНЫ ПРОМОТОРА cbh1")
    print("-" * 90)

    seq_len = len(sequence)
    print(f"\n1. Дистальный энхансерный регион (1-{seq_len//3}): Xyr1, ACE сайты")
    print(f"2. Проксимальный регуляторный регион ({seq_len//3}-{2*seq_len//3}): CRE1, HAP сайты")
    print(f"3. Коровый промотор ({2*seq_len//3}-{seq_len}): TATA-бокс, Inr, TSS")
    print(f"4. 5'-UTR ({seq_len}-ATG): Лидерная последовательность")

    # Специфичный анализ для cbh1
    print("\n" + "-" * 90)
    print("РЕГУЛЯЦИЯ ЭКСПРЕССИИ cbh1")
    print("-" * 90)

    # Подсчет ключевых элементов
    xyr1_count = sum(1 for e in results['elements'] if 'Xyr1' in e['name'])
    ace_count = sum(1 for e in results['elements'] if 'ACE' in e['name'])
    cre1_count = sum(1 for e in results['elements'] if 'CRE1' in e['name'])
    hap_count = sum(1 for e in results['elements'] if 'HAP' in e['name'])

    print(f"""
АКТИВАЦИЯ (индукция целлюлозой/софорозой):
• Xyr1 сайты: {xyr1_count} - главный активатор целлюлазных генов
• ACE сайты: {ace_count} - дополнительные активаторы
• HAP2/3/5 сайты: {hap_count} - CCAAT-связывающий комплекс

РЕПРЕССИЯ (глюкоза, CCR):
• CRE1 сайты: {cre1_count} - углеродная катаболитная репрессия

МЕХАНИЗМ ИНДУКЦИИ:
1. При отсутствии глюкозы и наличии целлюлозы
2. Xyr1 связывается с GGC(T/A)AA мотивами
3. Активируется транскрипция cbh1
4. При добавлении глюкозы CRE1 репрессирует промотор
""")

    return element_types


def create_visualization(sequence, results):
    """Создание визуализации промотора cbh1"""
    fig, axes = plt.subplots(4, 1, figsize=(18, 16))
    fig.suptitle('Анализ промотора cbh1 (Cel7A) из Trichoderma reesei', fontsize=14, fontweight='bold')

    # Специализированная цветовая схема для cbh1
    color_scheme = {
        'Xyr1': '#FF4444',      # Красный - главный активатор
        'ACE1': '#FF8C00',      # Оранжевый - активаторы
        'ACE2': '#FFD700',      # Золотой
        'ACE3': '#FFA500',      # Светло-оранжевый
        'CRE1': '#4169E1',      # Синий - репрессор
        'HAP': '#32CD32',       # Зеленый - HAP комплекс
        'TATA-box': '#9932CC',  # Фиолетовый
        'CAAT-box': '#4ECDC4',  # Бирюзовый
        'GC-box': '#45B7D1',    # Голубой
        'Initiator': '#FFEAA7', # Желтый
        'GATA': '#FFB6C1',      # Розовый
        'PacC': '#98D8C8',      # Мятный
        'CT-rich': '#96CEB4',   # Светло-зеленый
        'STRE': '#87CEEB',      # Небесно-голубой
        'Direct repeat': '#D3D3D3',  # Серый
        'Sophorose': '#FF69B4', # Ярко-розовый
        'Cellobiose': '#DA70D6', # Орхидея
    }

    seq_len = len(sequence)

    # === График 1: Карта регуляторных элементов ===
    ax1 = axes[0]
    ax1.set_xlim(0, seq_len + 50)
    ax1.set_ylim(-1, 10)

    # Рисуем основную линию ДНК
    ax1.plot([0, seq_len], [0, 0], 'k-', linewidth=3)
    ax1.plot([seq_len, seq_len + 30], [0, 0], 'g-', linewidth=3)

    # ATG
    ax1.annotate('ATG', xy=(seq_len + 15, 0), fontsize=10, ha='center', va='bottom',
                fontweight='bold', color='green')

    # Шкала
    for pos in range(0, seq_len, 200):
        ax1.plot([pos, pos], [-0.2, 0.2], 'k-', linewidth=1)
        ax1.text(pos, -0.5, str(pos+1), ha='center', fontsize=8)

    # Рисуем элементы
    y_levels = {}
    current_level = 1

    for elem in results['elements']:
        elem_type = elem['name'].split('(')[0].strip().split(' ')[0]
        color = color_scheme.get(elem_type, '#808080')

        # Определяем уровень y
        level = current_level
        for used_elem, used_level in y_levels.items():
            if used_level == level:
                if not (elem['end'] < used_elem[0] or elem['start'] > used_elem[1]):
                    level = max(y_levels.values()) + 1
                    if level > 9:
                        level = 1

        y_levels[(elem['start'], elem['end'])] = level

        # Прямоугольник элемента
        width = max(elem['end'] - elem['start'], 5)
        rect = FancyBboxPatch((elem['start'], level - 0.3), width, 0.6,
                              boxstyle="round,pad=0.02", facecolor=color,
                              edgecolor='black', linewidth=0.5, alpha=0.8)
        ax1.add_patch(rect)

        if width > 30:
            ax1.text(elem['start'] + width/2, level, elem_type[:8],
                    fontsize=6, ha='center', va='center')

    ax1.set_ylabel('Регуляторные элементы')
    ax1.set_title('Карта регуляторных элементов промотора cbh1')
    ax1.set_yticks([])
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['left'].set_visible(False)

    # Легенда (ключевые элементы для cbh1)
    key_elements = ['Xyr1', 'CRE1', 'ACE1', 'ACE2', 'HAP', 'TATA-box', 'CAAT-box', 'GC-box']
    legend_patches = [mpatches.Patch(color=color_scheme.get(name, '#808080'), label=name)
                      for name in key_elements]
    ax1.legend(handles=legend_patches, loc='upper right', fontsize=7, ncol=2)

    # === График 2: GC-содержание ===
    ax2 = axes[1]

    window_size = 50
    gc_profile = []
    positions = []

    for i in range(0, seq_len - window_size, 10):
        window = sequence[i:i+window_size]
        gc = (window.count('G') + window.count('C')) / window_size * 100
        gc_profile.append(gc)
        positions.append(i + window_size//2)

    ax2.fill_between(positions, gc_profile, alpha=0.3, color='blue')
    ax2.plot(positions, gc_profile, 'b-', linewidth=1.5)
    ax2.axhline(y=results['gc_content'], color='r', linestyle='--',
                label=f'Среднее GC: {results["gc_content"]:.1f}%')
    ax2.axhline(y=50, color='gray', linestyle=':', alpha=0.5)

    ax2.set_xlim(0, seq_len)
    ax2.set_ylim(0, 100)
    ax2.set_xlabel('Позиция (п.н.)')
    ax2.set_ylabel('GC-содержание (%)')
    ax2.set_title('Профиль GC-содержания (окно 50 п.н.)')
    ax2.legend(loc='upper right')
    ax2.grid(True, alpha=0.3)

    # === График 3: Плотность ключевых регуляторных элементов ===
    ax3 = axes[2]

    # Создаем отдельные профили для активаторов и репрессоров
    activator_density = np.zeros(seq_len)
    repressor_density = np.zeros(seq_len)

    for elem in results['elements']:
        for i in range(elem['start'], min(elem['end'], seq_len)):
            if any(act in elem['name'] for act in ['Xyr1', 'ACE1', 'ACE3', 'HAP', 'Sophorose']):
                activator_density[i] += 1
            elif any(rep in elem['name'] for rep in ['CRE1', 'ACE2']):
                repressor_density[i] += 1

    # Сглаживание
    try:
        from scipy.ndimage import gaussian_filter1d
        activator_smooth = gaussian_filter1d(activator_density, sigma=30)
        repressor_smooth = gaussian_filter1d(repressor_density, sigma=30)
    except ImportError:
        activator_smooth = np.convolve(activator_density, np.ones(60)/60, mode='same')
        repressor_smooth = np.convolve(repressor_density, np.ones(60)/60, mode='same')

    ax3.fill_between(range(seq_len), activator_smooth, alpha=0.4, color='green', label='Активаторы (Xyr1, ACE1, HAP)')
    ax3.fill_between(range(seq_len), -repressor_smooth, alpha=0.4, color='red', label='Репрессоры (CRE1, ACE2)')
    ax3.axhline(y=0, color='black', linewidth=0.5)

    ax3.set_xlim(0, seq_len)
    ax3.set_xlabel('Позиция (п.н.)')
    ax3.set_ylabel('Плотность элементов')
    ax3.set_title('Баланс активаторных и репрессорных элементов')
    ax3.legend(loc='upper right')
    ax3.grid(True, alpha=0.3)

    # === График 4: Схематичная структура промотора cbh1 ===
    ax4 = axes[3]
    ax4.set_xlim(-50, seq_len + 150)
    ax4.set_ylim(-3, 4)

    # Функциональные регионы
    regions = [
        (0, seq_len//3, 'Дистальный\nэнхансер\n(Xyr1, ACE)', '#FFCCCC'),
        (seq_len//3, 2*seq_len//3, 'Проксимальный\nрегион\n(CRE1, HAP)', '#CCE5FF'),
        (2*seq_len//3, seq_len, 'Коровый\nпромотор\n(TATA, Inr)', '#E0FFE0'),
        (seq_len, seq_len+40, "5'-UTR\n+ ATG", '#FFFACD'),
    ]

    for start, end, label, color in regions:
        rect = FancyBboxPatch((start, -0.7), end-start, 1.4,
                              boxstyle="round,pad=0.02", facecolor=color,
                              edgecolor='black', linewidth=1)
        ax4.add_patch(rect)
        ax4.text((start+end)/2, 0, label, ha='center', va='center', fontsize=8)

    # Стрелка транскрипции
    ax4.annotate('', xy=(seq_len + 120, 0), xytext=(seq_len + 60, 0),
                arrowprops=dict(arrowstyle='->', lw=2, color='green'))
    ax4.text(seq_len + 90, 0.8, 'Транскрипция', ha='center', fontsize=9, color='green')

    # Механизм регуляции
    ax4.text(seq_len//6, 2.5, '[+] Целлюлоза\nСофороза\nЛактоза', ha='center', fontsize=8,
             color='darkgreen', bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.5))
    ax4.text(seq_len//2, 2.5, '[-] Глюкоза\n(CCR)', ha='center', fontsize=8,
             color='darkred', bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.5))
    ax4.text(5*seq_len//6, 2.5, 'TSS\n(+1)', ha='center', fontsize=8,
             bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5))

    # Шкала позиций
    for pos in [0, -500, -1000, -1500]:
        abs_pos = seq_len + pos
        if abs_pos >= 0:
            ax4.plot([abs_pos, abs_pos], [-2.0, -1.7], 'k-', linewidth=1)
            ax4.text(abs_pos, -2.3, f'{pos}', ha='center', fontsize=8)

    ax4.text(seq_len//2, -2.8, 'Позиция относительно ATG (+1)', ha='center', fontsize=9)

    ax4.set_title('Схематичная структура промотора cbh1 и механизм регуляции')
    ax4.axis('off')

    plt.tight_layout()
    plt.savefig('pcbh1_promoter_analysis.png', dpi=150, bbox_inches='tight')
    plt.savefig('pcbh1_promoter_analysis.svg', format='svg', bbox_inches='tight')
    print("\n✓ Визуализация сохранена: pcbh1_promoter_analysis.png и pcbh1_promoter_analysis.svg")

    return fig


def create_sequence_annotation(sequence, results):
    """Создание аннотированной последовательности"""
    print("\n" + "=" * 90)
    print("АННОТИРОВАННАЯ ПОСЛЕДОВАТЕЛЬНОСТЬ")
    print("=" * 90)

    annotations = [' '] * len(sequence)

    # Маркеры для cbh1-специфичных элементов
    markers = {
        'Xyr1': 'X',
        'ACE1': '1',
        'ACE2': '2',
        'ACE3': '3',
        'CRE1': 'R',
        'HAP': 'H',
        'TATA': 'T',
        'CAAT': 'C',
        'GC-box': 'G',
        'Initiator': 'I',
        'GATA': 'A',
        'CT-rich': 'c',
    }

    for elem in results['elements']:
        for marker_name, marker_char in markers.items():
            if marker_name.lower() in elem['name'].lower():
                for i in range(elem['start'], min(elem['end'], len(sequence))):
                    annotations[i] = marker_char
                break

    # Вывод блоками по 60 символов
    block_size = 60
    print("\nЛегенда: X=Xyr1, 1=ACE1, 2=ACE2, 3=ACE3, R=CRE1, H=HAP, T=TATA, C=CAAT, G=GC, I=Inr, A=GATA")
    print("-" * 90)

    for i in range(0, len(sequence), block_size):
        block_seq = sequence[i:i+block_size]
        block_ann = ''.join(annotations[i:i+block_size])
        pos_from_atg = i - len(sequence)

        print(f"\n{i+1:4d}-{min(i+block_size, len(sequence)):4d} ({pos_from_atg:+5d})")
        print(f"     {block_seq}")
        print(f"     {block_ann}")

    print(f"\n{len(sequence)+1:4d}     ATG  (стартовый кодон)")


def main():
    """Главная функция"""
    sequence = PROMOTER_SEQUENCE

    print("\n" + "=" * 90)
    print("ПОСЛЕДОВАТЕЛЬНОСТЬ ПРОМОТОРА cbh1 (Trichoderma reesei)")
    print("=" * 90)
    print(f"\nДлина: {len(sequence)} п.н.")
    print("\n5'- " + sequence[:60] + "...")
    print("..." + sequence[-60:] + " -3' + ATG")

    # Анализ промотора
    results = analyze_promoter(sequence)

    # Вывод результатов
    element_types = print_analysis(results, sequence)

    # Аннотированная последовательность
    create_sequence_annotation(sequence, results)

    # Визуализация
    print("\n" + "=" * 90)
    print("СОЗДАНИЕ ВИЗУАЛИЗАЦИИ")
    print("=" * 90)

    try:
        fig = create_visualization(sequence, results)
        plt.close()
    except Exception as e:
        print(f"Ошибка при создании визуализации: {e}")
        import traceback
        traceback.print_exc()
        print("Попробуйте установить: pip install matplotlib scipy numpy")

    # Статистика
    print("\n" + "=" * 90)
    print("СВОДНАЯ СТАТИСТИКА")
    print("=" * 90)

    print(f"\nВсего идентифицировано элементов: {len(results['elements'])}")
    print("\nПо типам:")
    for elem_type, elements in sorted(element_types.items()):
        print(f"  • {elem_type}: {len(elements)}")

    # Сравнение с pGAP
    print("\n" + "=" * 90)
    print("СРАВНЕНИЕ cbh1 vs pGAP")
    print("=" * 90)
    print("""
┌──────────────────────┬────────────────────┬────────────────────┐
│ Характеристика       │ cbh1 (T. reesei)   │ pGAP (A. niger)    │
├──────────────────────┼────────────────────┼────────────────────┤
│ Тип промотора        │ Индуцибельный      │ Конститутивный     │
│ Индуктор             │ Целлюлоза, лактоза │ Глюкоза            │
│ Репрессор            │ Глюкоза (CCR)      │ -                  │
│ Уровень экспрессии   │ До 20-30 г/л       │ До 1-5 г/л         │
│ Ключевые TF          │ Xyr1, ACE, CRE1    │ Gcr1, CreA         │
│ Применение           │ Высокоуровневая    │ Конститутивная     │
│                      │ продукция белков   │ экспрессия         │
└──────────────────────┴────────────────────┴────────────────────┘
""")

    print("\n" + "=" * 90)
    print("ЗАКЛЮЧЕНИЕ")
    print("=" * 90)
    print("""
Промотор cbh1 (Cel7A) из T. reesei демонстрирует типичную структуру
индуцибельного промотора целлюлазных генов:

1. ИНДУКЦИЯ ЦЕЛЛЮЛОЗОЙ:
   - Множественные сайты Xyr1 (GGC[T/A]AA) - главный активатор
   - ACE1 сайты - дополнительная активация
   - Софорозо-респонсивные элементы
   - HAP2/3/5 комплекс (CCAAT) - общее усиление

2. РЕПРЕССИЯ ГЛЮКОЗОЙ (CCR):
   - CRE1 сайты (SYGGRG) - углеродная катаболитная репрессия
   - ACE2 сайты - дополнительная репрессия
   - При наличии глюкозы транскрипция блокируется

3. КОРОВЫЙ ПРОМОТОР:
   - TATA-бокс для позиционирования RNAPII
   - Инициаторные элементы (Inr)
   - Точка старта транскрипции (TSS)

БИОТЕХНОЛОГИЧЕСКОЕ ПРИМЕНЕНИЕ:
• Один из самых сильных известных промоторов (до 30 г/л белка)
• Используется для продукции промышленных ферментов
• Применяется для экспрессии рекомбинантных белков
• Требует индукции (лактоза - дешевый индуктор)
• Несовместим с глюкозой в среде (CCR)

РЕКОМЕНДАЦИИ ДЛЯ ИСПОЛЬЗОВАНИЯ:
1. Использовать среды без глюкозы
2. Индукция лактозой или целлюлозой
3. Двухфазная ферментация (рост → продукция)
4. Можно комбинировать с делецией CRE1 для снятия репрессии
""")


if __name__ == "__main__":
    main()
