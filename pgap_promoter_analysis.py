#!/usr/bin/env python3
"""
Анализ промотора pGAP (glyceraldehyde-3-phosphate dehydrogenase) из Aspergillus niger
Идентификация регуляторных элементов и визуализация структуры промотора
"""

import re
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, Arrow
import numpy as np

# Последовательность промотора pGAP из A. niger
PROMOTER_SEQUENCE = """
TACTACTATGAAAGACCGCGATGGGCCGATAGTATAGTTAGTTACTTCCATTACATCATCTCATCCGCCCGGTTCCTCGCCTCCGCGGCAGTCTACGGGTAGGATCGTAGCAAAAACCCGGGGGATAGACCCGTCGTCCCGAGCTGGAGTTCCGTATAACCTAGGTAGAAGGTATCAATTGAACCCGAACAACTGGCAAAACATTCTCGAGATCGTAGGAGTGAGTACCCGGCGTGATGGAGGGGGGAGCACGCTCATTGGTCCGTACGGCAGCTGCCGAGGGGGAGCAGGAGATCCAAATATCGTGAGTCTCCTGCTTTGCCCGGTGTATGAAACCGGAAAGGACTGCTGGGGAACTGGGGAGCGGCGCAAGCCGGGAATCCCAGCTGACAATTGACCCATCCTCATGCCGTGGCAGAGCTTGAGGTAGCTTTTGCCCCGTCTGTCTCCCCGGTGTGCGCATTCGACTGGGCGCGGCATCTGTGCCTCCTCCAGGAGCGGAGGACCCAGTAGTAAGTAGGCCTGACCTGGTCGTTGCGTCAGTCCAGAGGTTCCCTCCCCTACCCTTTTCTACTTCCCCTCCCCCGCCGCTCAACTTTTTCTTTTCCTTTTACTTTCTCTCTCTCTTCCTCTTCATCCATCCTCTCTTCATCACTTCCCTCTTCCCTTCATCCAATTCATCTTCCAAGTGAGTCTTCCTCCCCATCTGTCCCTCCATCTTTCCCATCATCATCTCCCCTCCCAGCTCCTCCCCTCCTCTCGTCTCCTCACGAAGCTTGACTAACCATTACCCCGCCACATAGACACATCTAAACA
""".replace('\n', '').replace(' ', '').upper()

# Стартовый кодон ATG (добавлен пользователем)
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
            'position_from_atg': pos - len(sequence)  # Позиция относительно ATG
        })
    return results

def analyze_promoter(sequence):
    """Комплексный анализ промотора"""
    results = {
        'length': len(sequence),
        'gc_content': (sequence.count('G') + sequence.count('C')) / len(sequence) * 100,
        'elements': []
    }

    # 1. TATA-бокс (консенсус: TATA(A/T)A(A/T), обычно -25 до -30 от TSS)
    # Расширенный паттерн для грибов
    tata_patterns = [
        (r'TATA[AT]A[AT]', 'TATA-box (canonical)'),
        (r'TATA[ATGC]A', 'TATA-box (variant)'),
        (r'TATAAA', 'TATA-box (TATAAA)'),
        (r'TATATA', 'TATA-box (TATATA)'),
    ]
    for pattern, name in tata_patterns:
        results['elements'].extend(find_pattern(sequence, pattern, name))

    # 2. CAAT-бокс (консенсус: CCAAT, обычно -70 до -80 от TSS)
    caat_patterns = [
        (r'CCAAT', 'CAAT-box'),
        (r'ATTGG', 'CAAT-box (reverse)'),
        (r'GG?CCAATC?', 'CAAT-box (extended)'),
    ]
    for pattern, name in caat_patterns:
        results['elements'].extend(find_pattern(sequence, pattern, name))

    # 3. GC-бокс (консенсус: GGGCGG или CCGCCC) - сайт связывания Sp1
    gc_patterns = [
        (r'GGGCGG', 'GC-box'),
        (r'CCGCCC', 'GC-box (reverse)'),
        (r'GG?GCGG?G', 'GC-box (variant)'),
    ]
    for pattern, name in gc_patterns:
        results['elements'].extend(find_pattern(sequence, pattern, name))

    # 4. CT-богатые регионы (характерны для грибных промоторов)
    ct_pattern = r'[CT]{8,}'
    results['elements'].extend(find_pattern(sequence, ct_pattern, 'CT-rich region'))

    # 5. Сайты инициатора (Inr) - PyPyAN(T/A)PyPy
    inr_patterns = [
        (r'[CT][CT]A[ATGC][TC][TC]', 'Initiator (Inr)'),
        (r'TCA[GC]T', 'Initiator (variant)'),
    ]
    for pattern, name in inr_patterns:
        results['elements'].extend(find_pattern(sequence, pattern, name))

    # 6. Сайты связывания транскрипционных факторов грибов

    # Gcr1 (glucose-responsive) - важен для GAP генов
    gcr1_patterns = [
        (r'CTTCC', 'Gcr1 binding site'),
        (r'GGAAG', 'Gcr1 binding site (reverse)'),
    ]
    for pattern, name in gcr1_patterns:
        results['elements'].extend(find_pattern(sequence, pattern, name))

    # Rap1 binding site (yeast-like, may be present in A. niger)
    rap1_pattern = r'[AT]CACCC[AT]'
    results['elements'].extend(find_pattern(sequence, rap1_pattern, 'Rap1-like site'))

    # CreA/Mig1 (carbon catabolite repression) - SYGGRG
    crea_patterns = [
        (r'[CG][CT]GG[AG]G', 'CreA/Mig1 binding site'),
        (r'GCGGGG', 'CreA binding site'),
    ]
    for pattern, name in crea_patterns:
        results['elements'].extend(find_pattern(sequence, pattern, name))

    # PacC (pH-responsive) - GCCARG
    pacc_pattern = r'GCCA[AG]G'
    results['elements'].extend(find_pattern(sequence, pacc_pattern, 'PacC binding site'))

    # AreA (nitrogen regulation) - GATA
    area_pattern = r'[AT]GATA[AG]'
    results['elements'].extend(find_pattern(sequence, area_pattern, 'AreA/GATA site'))

    # Hap complex binding (CCAAT-binding complex)
    hap_pattern = r'[CT]CAAT[CT]'
    results['elements'].extend(find_pattern(sequence, hap_pattern, 'Hap complex site'))

    # 7. Stress response elements
    stre_pattern = r'[AC]GGGG'  # STRE (stress response element)
    results['elements'].extend(find_pattern(sequence, stre_pattern, 'STRE (stress response)'))

    # 8. Heat shock element
    hse_pattern = r'[ATGC]GAA[ATGC]'
    results['elements'].extend(find_pattern(sequence, hse_pattern, 'HSE-like'))

    # 9. Kozak-подобная последовательность (около ATG)
    # В грибах: (A/C)A(A/C)ATGG
    # Проверяем конец последовательности
    if sequence[-10:]:
        kozak_match = re.search(r'[AC]A[AC]A[CT][AC]ATG', sequence[-20:] + 'ATG', re.IGNORECASE)
        if kozak_match:
            results['elements'].append({
                'name': 'Kozak-like sequence',
                'pattern': '[AC]A[AC]A[CT][AC]ATG',
                'match': kozak_match.group(),
                'start': len(sequence) - 20 + kozak_match.start(),
                'end': len(sequence),
                'position_from_atg': -20 + kozak_match.start()
            })

    # 10. Поиск повторяющихся элементов
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
            # Поиск повторов этого мотива
            count = sequence.count(motif)
            if count >= 2 and len(set(motif)) > 2:  # Минимум 2 повтора, не гомополимер
                positions = [m.start() for m in re.finditer(re.escape(motif), sequence)]
                if len(positions) >= 2 and positions[0] == i:
                    # Проверяем, не перекрывается ли с уже найденным
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
    return repeats[:10]  # Возвращаем только топ-10

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
    print("=" * 80)
    print("АНАЛИЗ ПРОМОТОРА pGAP ИЗ Aspergillus niger")
    print("=" * 80)
    print(f"\nДлина последовательности: {results['length']} п.н. (+ ATG)")
    print(f"GC-содержание: {results['gc_content']:.1f}%")
    print(f"AT-содержание: {100 - results['gc_content']:.1f}%")

    print("\n" + "-" * 80)
    print("ИДЕНТИФИЦИРОВАННЫЕ РЕГУЛЯТОРНЫЕ ЭЛЕМЕНТЫ")
    print("-" * 80)

    # Группировка элементов по типу
    element_types = {}
    for elem in results['elements']:
        elem_type = elem['name'].split('(')[0].strip()
        if elem_type not in element_types:
            element_types[elem_type] = []
        element_types[elem_type].append(elem)

    for elem_type, elements in sorted(element_types.items()):
        print(f"\n{elem_type}:")
        for elem in elements:
            pos_str = f"{elem['position_from_atg']:+d}" if elem['position_from_atg'] < 0 else f"+{elem['position_from_atg']}"
            print(f"  • Позиция {elem['start']+1}-{elem['end']} ({pos_str} от ATG): {elem['match']}")

    print("\n" + "-" * 80)
    print("ФУНКЦИОНАЛЬНЫЕ РЕГИОНЫ ПРОМОТОРА")
    print("-" * 80)

    # Определение функциональных регионов
    print(f"\n1. Дистальный регион (1-{len(sequence)//3}): Энхансерные элементы")
    print(f"2. Проксимальный регион ({len(sequence)//3}-{2*len(sequence)//3}): Регуляторные элементы")
    print(f"3. Коровый регион ({2*len(sequence)//3}-{len(sequence)}): TATA-бокс, Inr, TSS")
    print(f"4. 5'-UTR ({len(sequence)}-ATG): Лидерная последовательность")

    # Анализ корового региона
    core_region = sequence[-150:]
    print("\n" + "-" * 80)
    print("АНАЛИЗ КОРОВОГО РЕГИОНА (последние 150 п.н.)")
    print("-" * 80)

    # Поиск потенциального TSS
    pyrimidine_rich = []
    for i in range(len(core_region) - 10):
        window = core_region[i:i+10]
        py_count = window.count('C') + window.count('T')
        if py_count >= 7:
            pyrimidine_rich.append((i + len(sequence) - 150, py_count))

    if pyrimidine_rich:
        print("\nПиримидин-богатые регионы (потенциальные TSS):")
        for pos, count in pyrimidine_rich[:5]:
            print(f"  • Позиция {pos+1}: {count}/10 пиримидинов")

    return element_types

def create_visualization(sequence, results):
    """Создание визуализации промотора"""
    fig, axes = plt.subplots(4, 1, figsize=(16, 14))
    fig.suptitle('Анализ промотора pGAP из Aspergillus niger', fontsize=14, fontweight='bold')

    # Цветовая схема для элементов
    color_scheme = {
        'TATA-box': '#FF6B6B',
        'CAAT-box': '#4ECDC4',
        'GC-box': '#45B7D1',
        'CT-rich region': '#96CEB4',
        'Initiator': '#FFEAA7',
        'Gcr1': '#DDA0DD',
        'CreA': '#F0E68C',
        'PacC': '#98D8C8',
        'AreA': '#FFB6C1',
        'STRE': '#87CEEB',
        'Hap': '#DEB887',
        'Direct repeat': '#D3D3D3',
        'HSE': '#FFA07A',
        'Kozak': '#90EE90',
        'Rap1': '#E6E6FA',
    }

    seq_len = len(sequence)

    # === График 1: Карта регуляторных элементов ===
    ax1 = axes[0]
    ax1.set_xlim(0, seq_len + 50)
    ax1.set_ylim(-1, 8)

    # Рисуем основную линию ДНК
    ax1.plot([0, seq_len], [0, 0], 'k-', linewidth=3)
    ax1.plot([seq_len, seq_len + 30], [0, 0], 'g-', linewidth=3)  # ATG

    # Добавляем ATG
    ax1.annotate('ATG', xy=(seq_len + 15, 0), fontsize=10, ha='center', va='bottom',
                fontweight='bold', color='green')

    # Добавляем шкалу
    for pos in range(0, seq_len, 100):
        ax1.plot([pos, pos], [-0.2, 0.2], 'k-', linewidth=1)
        ax1.text(pos, -0.5, str(pos+1), ha='center', fontsize=8)

    # Рисуем элементы
    y_levels = {}
    current_level = 1

    for elem in results['elements']:
        elem_type = elem['name'].split('(')[0].strip().split(' ')[0]
        color = color_scheme.get(elem_type, '#808080')

        # Определяем уровень y для избежания перекрытий
        level = current_level
        for used_elem, used_level in y_levels.items():
            if used_level == level:
                if not (elem['end'] < used_elem[0] or elem['start'] > used_elem[1]):
                    level = max(y_levels.values()) + 1

        y_levels[(elem['start'], elem['end'])] = level

        # Рисуем прямоугольник элемента
        width = max(elem['end'] - elem['start'], 3)
        rect = FancyBboxPatch((elem['start'], level - 0.3), width, 0.6,
                              boxstyle="round,pad=0.02", facecolor=color,
                              edgecolor='black', linewidth=0.5, alpha=0.8)
        ax1.add_patch(rect)

        # Добавляем метку
        if width > 20:
            ax1.text(elem['start'] + width/2, level, elem_type[:10],
                    fontsize=6, ha='center', va='center')

    ax1.set_ylabel('Регуляторные элементы')
    ax1.set_title('Карта регуляторных элементов промотора')
    ax1.set_yticks([])
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['left'].set_visible(False)

    # Легенда
    legend_patches = [mpatches.Patch(color=color, label=name)
                      for name, color in list(color_scheme.items())[:8]]
    ax1.legend(handles=legend_patches, loc='upper right', fontsize=7, ncol=2)

    # === График 2: GC-содержание по длине ===
    ax2 = axes[1]

    window_size = 30
    gc_profile = []
    positions = []

    for i in range(0, seq_len - window_size, 5):
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
    ax2.set_title('Профиль GC-содержания (окно 30 п.н.)')
    ax2.legend(loc='upper right')
    ax2.grid(True, alpha=0.3)

    # === График 3: Плотность элементов ===
    ax3 = axes[2]

    # Подсчет элементов в окнах
    element_density = np.zeros(seq_len)
    for elem in results['elements']:
        for i in range(elem['start'], min(elem['end'], seq_len)):
            element_density[i] += 1

    # Сглаживание
    from scipy.ndimage import gaussian_filter1d
    try:
        smoothed_density = gaussian_filter1d(element_density, sigma=20)
    except:
        smoothed_density = np.convolve(element_density, np.ones(40)/40, mode='same')

    ax3.fill_between(range(seq_len), smoothed_density, alpha=0.4, color='purple')
    ax3.plot(range(seq_len), smoothed_density, 'purple', linewidth=1.5)

    ax3.set_xlim(0, seq_len)
    ax3.set_xlabel('Позиция (п.н.)')
    ax3.set_ylabel('Плотность элементов')
    ax3.set_title('Плотность регуляторных элементов')
    ax3.grid(True, alpha=0.3)

    # === График 4: Схематичная структура ===
    ax4 = axes[3]
    ax4.set_xlim(-50, seq_len + 100)
    ax4.set_ylim(-2, 3)

    # Функциональные регионы
    regions = [
        (0, seq_len//3, 'Дистальный\nрегион', '#FFE4E1'),
        (seq_len//3, 2*seq_len//3, 'Проксимальный\nрегион', '#E6E6FA'),
        (2*seq_len//3, seq_len, 'Коровый\nрегион', '#E0FFE0'),
        (seq_len, seq_len+30, "5'-UTR\n+ ATG", '#FFFACD'),
    ]

    for start, end, label, color in regions:
        rect = FancyBboxPatch((start, -0.5), end-start, 1,
                              boxstyle="round,pad=0.02", facecolor=color,
                              edgecolor='black', linewidth=1)
        ax4.add_patch(rect)
        ax4.text((start+end)/2, 0, label, ha='center', va='center', fontsize=9)

    # Стрелка направления транскрипции
    ax4.annotate('', xy=(seq_len + 80, 0), xytext=(seq_len + 40, 0),
                arrowprops=dict(arrowstyle='->', lw=2, color='green'))
    ax4.text(seq_len + 60, 0.8, 'Транскрипция', ha='center', fontsize=9, color='green')

    # Добавляем ключевые элементы
    ax4.text(seq_len//6, 1.5, 'Энхансеры\nTF сайты', ha='center', fontsize=8, style='italic')
    ax4.text(seq_len//2, 1.5, 'CAAT-бокс\nGC-бокс\nCreA, PacC', ha='center', fontsize=8, style='italic')
    ax4.text(5*seq_len//6, 1.5, 'TATA-бокс\nInr\nTSS', ha='center', fontsize=8, style='italic')

    # Шкала позиций относительно ATG
    for pos in [0, -200, -400, -600, -800]:
        abs_pos = seq_len + pos
        if abs_pos >= 0:
            ax4.plot([abs_pos, abs_pos], [-1.5, -1.2], 'k-', linewidth=1)
            ax4.text(abs_pos, -1.8, f'{pos}', ha='center', fontsize=8)

    ax4.text(seq_len//2, -2.3, 'Позиция относительно ATG (+1)', ha='center', fontsize=9)

    ax4.set_title('Схематичная структура промотора pGAP')
    ax4.axis('off')

    plt.tight_layout()
    plt.savefig('pgap_promoter_analysis.png', dpi=150, bbox_inches='tight')
    plt.savefig('pgap_promoter_analysis.svg', format='svg', bbox_inches='tight')
    print("\n✓ Визуализация сохранена: pgap_promoter_analysis.png и pgap_promoter_analysis.svg")

    return fig

def create_sequence_annotation(sequence, results):
    """Создание аннотированной последовательности"""
    print("\n" + "=" * 80)
    print("АННОТИРОВАННАЯ ПОСЛЕДОВАТЕЛЬНОСТЬ")
    print("=" * 80)

    # Создаем массив аннотаций
    annotations = [' '] * len(sequence)

    # Маркируем элементы
    markers = {
        'TATA': 'T',
        'CAAT': 'C',
        'GC-box': 'G',
        'CT-rich': 'R',
        'Initiator': 'I',
        'Gcr1': '1',
        'CreA': 'A',
        'PacC': 'P',
        'STRE': 'S',
    }

    for elem in results['elements']:
        for marker_name, marker_char in markers.items():
            if marker_name.lower() in elem['name'].lower():
                for i in range(elem['start'], min(elem['end'], len(sequence))):
                    annotations[i] = marker_char
                break

    # Выводим последовательность блоками по 60 символов
    block_size = 60
    print("\nЛегенда: T=TATA, C=CAAT, G=GC-box, R=CT-rich, I=Inr, 1=Gcr1, A=CreA, P=PacC, S=STRE")
    print("-" * 80)

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

    print("\n" + "=" * 80)
    print("ПОСЛЕДОВАТЕЛЬНОСТЬ ПРОМОТОРА pGAP (Aspergillus niger)")
    print("=" * 80)
    print(f"\nДлина: {len(sequence)} п.н.")
    print("\n5'- " + sequence[:50] + "...")
    print("..." + sequence[-50:] + " -3' + ATG")

    # Анализ промотора
    results = analyze_promoter(sequence)

    # Вывод результатов
    element_types = print_analysis(results, sequence)

    # Аннотированная последовательность
    create_sequence_annotation(sequence, results)

    # Визуализация
    print("\n" + "=" * 80)
    print("СОЗДАНИЕ ВИЗУАЛИЗАЦИИ")
    print("=" * 80)

    try:
        fig = create_visualization(sequence, results)
        plt.close()
    except Exception as e:
        print(f"Ошибка при создании визуализации: {e}")
        print("Попробуйте установить: pip install matplotlib scipy numpy")

    # Статистика
    print("\n" + "=" * 80)
    print("СВОДНАЯ СТАТИСТИКА")
    print("=" * 80)

    print(f"\nВсего идентифицировано элементов: {len(results['elements'])}")
    print("\nПо типам:")
    for elem_type, elements in sorted(element_types.items()):
        print(f"  • {elem_type}: {len(elements)}")

    print("\n" + "=" * 80)
    print("ЗАКЛЮЧЕНИЕ")
    print("=" * 80)
    print("""
Промотор pGAP из A. niger демонстрирует типичную структуру грибного промотора:

1. КОРОВЫЙ ПРОМОТОР:
   - TATA-подобные элементы для позиционирования РНК-полимеразы II
   - Инициаторные последовательности (Inr) около TSS
   - Пиримидин-богатые регионы

2. РЕГУЛЯТОРНЫЕ ЭЛЕМЕНТЫ:
   - Gcr1 сайты - характерны для гликолитических генов
   - CreA/Mig1 сайты - углеродная катаболитная репрессия
   - CAAT-боксы - общие усилители транскрипции

3. ОСОБЕННОСТИ:
   - Умеренное GC-содержание (типично для A. niger)
   - CT-богатые регионы - могут влиять на стабильность мРНК
   - Множественные регуляторные сайты - обеспечивают тонкую регуляцию

Этот промотор широко используется для гетерологичной экспрессии генов
в промышленных штаммах Aspergillus благодаря его конститутивной
и высокой активности при росте на глюкозе.
""")

if __name__ == "__main__":
    main()
