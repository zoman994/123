#!/usr/bin/env python3
"""
Анализ возможностей усиления промотора pGAP и снижения репрессии
Стратегии модификации для увеличения экспрессии генов

Автор: Анализ на основе данных промотора pGAP из Aspergillus niger
"""

import re
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, Rectangle, FancyArrowPatch
import numpy as np

# Последовательность промотора pGAP из A. niger
PROMOTER_SEQUENCE = """
TACTACTATGAAAGACCGCGATGGGCCGATAGTATAGTTAGTTACTTCCATTACATCATCTCATCCGCCCGGTTCCTCGCCTCCGCGGCAGTCTACGGGTAGGATCGTAGCAAAAACCCGGGGGATAGACCCGTCGTCCCGAGCTGGAGTTCCGTATAACCTAGGTAGAAGGTATCAATTGAACCCGAACAACTGGCAAAACATTCTCGAGATCGTAGGAGTGAGTACCCGGCGTGATGGAGGGGGGAGCACGCTCATTGGTCCGTACGGCAGCTGCCGAGGGGGAGCAGGAGATCCAAATATCGTGAGTCTCCTGCTTTGCCCGGTGTATGAAACCGGAAAGGACTGCTGGGGAACTGGGGAGCGGCGCAAGCCGGGAATCCCAGCTGACAATTGACCCATCCTCATGCCGTGGCAGAGCTTGAGGTAGCTTTTGCCCCGTCTGTCTCCCCGGTGTGCGCATTCGACTGGGCGCGGCATCTGTGCCTCCTCCAGGAGCGGAGGACCCAGTAGTAAGTAGGCCTGACCTGGTCGTTGCGTCAGTCCAGAGGTTCCCTCCCCTACCCTTTTCTACTTCCCCTCCCCCGCCGCTCAACTTTTTCTTTTCCTTTTACTTTCTCTCTCTCTTCCTCTTCATCCATCCTCTCTTCATCACTTCCCTCTTCCCTTCATCCAATTCATCTTCCAAGTGAGTCTTCCTCCCCATCTGTCCCTCCATCTTTCCCATCATCATCTCCCCTCCCAGCTCCTCCCCTCCTCTCGTCTCCTCACGAAGCTTGACTAACCATTACCCCGCCACATAGACACATCTAAACA
""".replace('\n', '').replace(' ', '').upper()


class PromoterEnhancementAnalyzer:
    """Класс для анализа возможностей усиления промотора"""

    def __init__(self, sequence):
        self.sequence = sequence
        self.length = len(sequence)
        self.elements = self._find_all_elements()
        self.repressor_sites = self._identify_repressor_sites()
        self.activator_sites = self._identify_activator_sites()

    def _find_pattern(self, pattern, name):
        """Поиск паттерна в последовательности"""
        results = []
        for match in re.finditer(pattern, self.sequence, re.IGNORECASE):
            results.append({
                'name': name,
                'pattern': pattern,
                'match': match.group(),
                'start': match.start(),
                'end': match.end(),
                'position_from_atg': match.start() - self.length
            })
        return results

    def _find_all_elements(self):
        """Поиск всех регуляторных элементов"""
        elements = []

        # Базовые элементы
        patterns = [
            (r'TATA[AT]A[AT]', 'TATA-box'),
            (r'TATATA', 'TATA-box'),
            (r'CCAAT', 'CAAT-box'),
            (r'ATTGG', 'CAAT-box (reverse)'),
            (r'GGGCGG', 'GC-box'),
            (r'CCGCCC', 'GC-box (reverse)'),
            (r'CTTCC', 'Gcr1 site'),
            (r'GGAAG', 'Gcr1 site (reverse)'),
            (r'[CG][CT]GG[AG]G', 'CreA/Mig1 site'),
            (r'GCGGGG', 'CreA site'),
            (r'GCCA[AG]G', 'PacC site'),
            (r'[AT]GATA[AG]', 'AreA/GATA site'),
            (r'[CT]CAAT[CT]', 'Hap complex site'),
            (r'[AC]GGGG', 'STRE'),
            (r'[CT]{8,}', 'CT-rich region'),
        ]

        for pattern, name in patterns:
            elements.extend(self._find_pattern(pattern, name))

        return elements

    def _identify_repressor_sites(self):
        """Идентификация сайтов репрессоров"""
        repressors = []

        # CreA/Mig1 - углеродная катаболитная репрессия
        crea_patterns = [
            (r'[CG][CT]GG[AG]G', 'CreA/Mig1', 'Carbon catabolite repression'),
            (r'GCGGGG', 'CreA', 'Carbon catabolite repression'),
            (r'SYGGRG', 'CreA consensus', 'Glucose repression'),
        ]

        for pattern, name, function in crea_patterns:
            for match in re.finditer(pattern, self.sequence, re.IGNORECASE):
                repressors.append({
                    'name': name,
                    'function': function,
                    'match': match.group(),
                    'start': match.start(),
                    'end': match.end(),
                    'position_from_atg': match.start() - self.length,
                    'modification_strategy': 'DELETE or MUTATE',
                    'expected_effect': 'Relief from carbon catabolite repression',
                    'risk_level': 'MEDIUM'
                })

        # Потенциальные сайты репрессии в AT-богатых регионах
        at_rich = re.finditer(r'[AT]{10,}', self.sequence)
        for match in at_rich:
            repressors.append({
                'name': 'AT-rich region',
                'function': 'Potential silencer/nucleosome positioning',
                'match': match.group()[:20] + '...' if len(match.group()) > 20 else match.group(),
                'start': match.start(),
                'end': match.end(),
                'position_from_atg': match.start() - self.length,
                'modification_strategy': 'REPLACE with neutral sequence',
                'expected_effect': 'May reduce nucleosome-mediated repression',
                'risk_level': 'LOW'
            })

        return repressors

    def _identify_activator_sites(self):
        """Идентификация сайтов активаторов"""
        activators = []

        # Gcr1 - глюкозо-респонсивная активация гликолитических генов
        gcr1_patterns = [
            (r'CTTCC', 'Gcr1', 'Glycolytic gene activation'),
            (r'GGAAG', 'Gcr1 (reverse)', 'Glycolytic gene activation'),
        ]

        for pattern, name, function in gcr1_patterns:
            for match in re.finditer(pattern, self.sequence, re.IGNORECASE):
                activators.append({
                    'name': name,
                    'function': function,
                    'match': match.group(),
                    'start': match.start(),
                    'end': match.end(),
                    'position_from_atg': match.start() - self.length,
                    'enhancement_strategy': 'DUPLICATE or optimize spacing',
                    'expected_effect': '2-5x increase in expression'
                })

        # CAAT-box
        for match in re.finditer(r'CCAAT|ATTGG', self.sequence, re.IGNORECASE):
            activators.append({
                'name': 'CAAT-box',
                'function': 'General transcription enhancement',
                'match': match.group(),
                'start': match.start(),
                'end': match.end(),
                'position_from_atg': match.start() - self.length,
                'enhancement_strategy': 'ADD consensus CAAT-box at -80',
                'expected_effect': '1.5-3x increase'
            })

        # GC-box (Sp1 binding)
        for match in re.finditer(r'GGGCGG|CCGCCC', self.sequence, re.IGNORECASE):
            activators.append({
                'name': 'GC-box',
                'function': 'Sp1-mediated activation',
                'match': match.group(),
                'start': match.start(),
                'end': match.end(),
                'position_from_atg': match.start() - self.length,
                'enhancement_strategy': 'ADD multiple GC-boxes',
                'expected_effect': '2-4x increase'
            })

        return activators

    def analyze_enhancement_strategies(self):
        """Анализ стратегий усиления промотора"""
        strategies = {
            'core_promoter': [],
            'upstream_activation': [],
            'repression_relief': [],
            'hybrid_promoter': [],
            'synthetic_elements': []
        }

        # 1. Оптимизация корового промотора
        strategies['core_promoter'] = [
            {
                'strategy': 'TATA-box optimization',
                'description': 'Замена существующих TATA-подобных элементов на консенсусную последовательность TATAAAA',
                'position': '-30 to -25 от TSS',
                'current_analysis': self._analyze_tata_region(),
                'modification': 'Мутировать в TATAAAA для максимальной эффективности инициации',
                'expected_effect': '2-5x увеличение базальной транскрипции',
                'risk': 'LOW - хорошо изученная модификация',
                'priority': 'HIGH'
            },
            {
                'strategy': 'Initiator (Inr) optimization',
                'description': 'Оптимизация последовательности инициатора транскрипции',
                'position': '-2 to +5 относительно TSS',
                'modification': 'Создать консенсус PyPyAN(T/A)PyPy',
                'expected_effect': '1.5-2x увеличение',
                'risk': 'LOW',
                'priority': 'MEDIUM'
            },
            {
                'strategy': 'Kozak sequence optimization',
                'description': 'Оптимизация контекста стартового кодона для улучшения инициации трансляции',
                'position': '-6 to +4 относительно ATG',
                'modification': 'Создать последовательность (A/G)CCATGG',
                'expected_effect': '2-10x увеличение эффективности трансляции',
                'risk': 'LOW',
                'priority': 'HIGH'
            }
        ]

        # 2. Upstream Activating Sequences (UAS)
        strategies['upstream_activation'] = [
            {
                'strategy': 'Gcr1 site multiplication',
                'description': 'Добавление дополнительных сайтов связывания Gcr1',
                'current_count': len([e for e in self.activator_sites if 'Gcr1' in e['name']]),
                'modification': 'Добавить 2-4 тандемных Gcr1 сайта (CTTCC) в оптимальном спейсинге',
                'optimal_position': '-200 to -400 от ATG',
                'expected_effect': '3-8x увеличение на глюкозе',
                'risk': 'MEDIUM',
                'priority': 'HIGH'
            },
            {
                'strategy': 'UAS element insertion',
                'description': 'Вставка проверенных UAS элементов из других промоторов',
                'options': [
                    'UAS из pAOX1 (для индуцибельности метанолом)',
                    'UAS из pGAL (для галактозной индукции)',
                    'UAS из pENO (конститутивная активация)'
                ],
                'modification': 'Клонировать UAS перед коровым промотором',
                'expected_effect': '5-50x увеличение (зависит от UAS)',
                'risk': 'MEDIUM-HIGH',
                'priority': 'HIGH'
            },
            {
                'strategy': 'Enhancer addition',
                'description': 'Добавление энхансерных последовательностей',
                'modification': 'Вставка CMV энхансера или грибного энхансера',
                'expected_effect': '10-100x увеличение',
                'risk': 'HIGH - может влиять на регуляцию',
                'priority': 'MEDIUM'
            }
        ]

        # 3. Снятие репрессии
        crea_sites = [r for r in self.repressor_sites if 'CreA' in r['name'] or 'Mig1' in r['name']]
        strategies['repression_relief'] = [
            {
                'strategy': 'CreA/Mig1 site deletion',
                'description': 'Удаление сайтов углеродной катаболитной репрессии',
                'current_sites': len(crea_sites),
                'sites_positions': [f"{s['position_from_atg']:+d}" for s in crea_sites],
                'modification': 'Делеция или мутация последовательностей SYGGRG',
                'mutation_options': [
                    'GCGGGG → GCAGAG (сохранение GC, потеря связывания)',
                    'CTGGAG → CTGAAG (минимальная мутация)',
                    'Полная делеция 6-8 п.н.'
                ],
                'expected_effect': '2-10x увеличение в присутствии глюкозы',
                'risk': 'MEDIUM - проверить экспрессию при разных источниках углерода',
                'priority': 'HIGH'
            },
            {
                'strategy': 'CreA host deletion',
                'description': 'Использование штамма с делецией гена creA',
                'modification': 'Нокаут creA в штамме-хозяине',
                'expected_effect': '5-20x увеличение для репрессированных генов',
                'risk': 'HIGH - плейотропные эффекты на метаболизм',
                'priority': 'MEDIUM'
            },
            {
                'strategy': 'Derepression by glucose limitation',
                'description': 'Культивирование в условиях лимитированной глюкозы',
                'modification': 'Fed-batch с контролем глюкозы <0.1%',
                'expected_effect': 'Полное снятие катаболитной репрессии',
                'risk': 'LOW - стандартная практика',
                'priority': 'LOW'
            }
        ]

        # 4. Гибридные промоторы
        strategies['hybrid_promoter'] = [
            {
                'strategy': 'pGAP-pAOX1 hybrid',
                'description': 'Комбинация UAS от pAOX1 с коровым pGAP',
                'modification': 'UAS(AOX1) + Core(GAP)',
                'expected_effect': 'Метанол-индуцибельный + высокая базальная активность',
                'applications': 'Двухфазная экспрессия',
                'risk': 'MEDIUM',
                'priority': 'MEDIUM'
            },
            {
                'strategy': 'Tandem promoter',
                'description': 'Тандемное расположение двух промоторов',
                'modification': 'pGAP-pGAP или pGAP-pTEF тандем',
                'expected_effect': '2-4x увеличение',
                'risk': 'LOW',
                'priority': 'MEDIUM'
            },
            {
                'strategy': 'Bidirectional promoter',
                'description': 'Создание бидирекционального промотора',
                'modification': 'Head-to-head конфигурация с pGAP',
                'expected_effect': 'Экспрессия двух генов с одного промотора',
                'risk': 'MEDIUM',
                'priority': 'LOW'
            }
        ]

        # 5. Синтетические элементы
        strategies['synthetic_elements'] = [
            {
                'strategy': 'Synthetic TFBS array',
                'description': 'Синтетический массив сайтов связывания ТФ',
                'modification': '''Создание синтетической кассеты:
                    [Gcr1]x4 - [CAAT] - [GC-box]x2 - TATA - Inr - ATG''',
                'expected_effect': '10-100x увеличение возможно',
                'risk': 'MEDIUM-HIGH',
                'priority': 'MEDIUM'
            },
            {
                'strategy': '5\'-UTR engineering',
                'description': 'Оптимизация 5\'-нетранслируемой области',
                'modifications': [
                    'Удаление вторичных структур (uORF)',
                    'Оптимизация длины (оптимум 50-100 нт)',
                    'Добавление IRES для cap-независимой трансляции'
                ],
                'expected_effect': '2-5x увеличение трансляции',
                'risk': 'LOW-MEDIUM',
                'priority': 'HIGH'
            },
            {
                'strategy': 'Codon optimization context',
                'description': 'Оптимизация первых кодонов после ATG',
                'modification': 'A-rich последовательность после ATG улучшает инициацию',
                'expected_effect': '1.5-3x увеличение',
                'risk': 'LOW',
                'priority': 'MEDIUM'
            }
        ]

        return strategies

    def _analyze_tata_region(self):
        """Анализ TATA-региона"""
        # TATA-бокс обычно находится в позиции -30 до -25 от TSS
        # Предполагаем TSS примерно -50 от ATG для грибов
        tata_region = self.sequence[-100:-50]  # Регион где ожидается TATA

        tata_elements = []
        for pattern in [r'TATA[AT]A', r'TATATA', r'TATAAA']:
            for match in re.finditer(pattern, tata_region, re.IGNORECASE):
                tata_elements.append({
                    'sequence': match.group(),
                    'position_in_region': match.start(),
                    'position_from_atg': match.start() - 100,
                    'quality': 'optimal' if match.group() in ['TATAAA', 'TATATA'] else 'suboptimal'
                })

        return {
            'found': len(tata_elements),
            'elements': tata_elements,
            'recommendation': 'Add consensus TATAAA if none found' if not tata_elements else 'Optimize existing'
        }

    def calculate_modification_priority(self):
        """Расчет приоритета модификаций"""
        priorities = []

        strategies = self.analyze_enhancement_strategies()

        for category, items in strategies.items():
            for item in items:
                if isinstance(item, dict) and 'priority' in item:
                    score = {'HIGH': 3, 'MEDIUM': 2, 'LOW': 1}.get(item['priority'], 0)
                    risk_penalty = {'HIGH': -1, 'MEDIUM-HIGH': -0.5, 'MEDIUM': 0, 'LOW-MEDIUM': 0.25, 'LOW': 0.5}
                    risk = item.get('risk', 'MEDIUM').split(' - ')[0]
                    risk_score = risk_penalty.get(risk, 0)

                    priorities.append({
                        'category': category,
                        'strategy': item['strategy'],
                        'priority_score': score + risk_score,
                        'expected_effect': item.get('expected_effect', 'Unknown'),
                        'risk': item.get('risk', 'Unknown')
                    })

        return sorted(priorities, key=lambda x: x['priority_score'], reverse=True)

    def generate_report(self):
        """Генерация полного отчета"""
        report = []
        report.append("=" * 90)
        report.append("АНАЛИЗ ВОЗМОЖНОСТЕЙ УСИЛЕНИЯ ПРОМОТОРА pGAP")
        report.append("=" * 90)

        report.append(f"\nДлина промотора: {self.length} п.н.")
        report.append(f"Всего регуляторных элементов: {len(self.elements)}")
        report.append(f"Сайтов репрессоров: {len(self.repressor_sites)}")
        report.append(f"Сайтов активаторов: {len(self.activator_sites)}")

        # Репрессорные сайты
        report.append("\n" + "-" * 90)
        report.append("ИДЕНТИФИЦИРОВАННЫЕ САЙТЫ РЕПРЕССИИ (мишени для модификации)")
        report.append("-" * 90)

        for rep in self.repressor_sites:
            report.append(f"\n  ► {rep['name']} [{rep['function']}]")
            report.append(f"    Позиция: {rep['position_from_atg']:+d} от ATG")
            report.append(f"    Последовательность: {rep['match']}")
            report.append(f"    Стратегия: {rep['modification_strategy']}")
            report.append(f"    Ожидаемый эффект: {rep['expected_effect']}")
            report.append(f"    Риск: {rep['risk_level']}")

        # Активаторные сайты
        report.append("\n" + "-" * 90)
        report.append("САЙТЫ АКТИВАЦИИ (можно усилить)")
        report.append("-" * 90)

        for act in self.activator_sites:
            report.append(f"\n  ► {act['name']} [{act['function']}]")
            report.append(f"    Позиция: {act['position_from_atg']:+d} от ATG")
            report.append(f"    Последовательность: {act['match']}")
            report.append(f"    Стратегия усиления: {act['enhancement_strategy']}")
            report.append(f"    Ожидаемый эффект: {act['expected_effect']}")

        # Стратегии модификации
        strategies = self.analyze_enhancement_strategies()

        for category, items in strategies.items():
            report.append("\n" + "=" * 90)
            report.append(f"СТРАТЕГИИ: {category.upper().replace('_', ' ')}")
            report.append("=" * 90)

            for item in items:
                if isinstance(item, dict):
                    report.append(f"\n  ▶ {item.get('strategy', 'Unknown')}")
                    report.append(f"    {item.get('description', '')}")
                    if 'modification' in item:
                        report.append(f"    Модификация: {item['modification']}")
                    if 'expected_effect' in item:
                        report.append(f"    Эффект: {item['expected_effect']}")
                    if 'risk' in item:
                        report.append(f"    Риск: {item['risk']}")
                    if 'priority' in item:
                        report.append(f"    Приоритет: {item['priority']}")

        # Приоритеты
        report.append("\n" + "=" * 90)
        report.append("РЕКОМЕНДОВАННЫЙ ПОРЯДОК МОДИФИКАЦИЙ")
        report.append("=" * 90)

        priorities = self.calculate_modification_priority()
        for i, p in enumerate(priorities[:10], 1):
            report.append(f"\n  {i}. {p['strategy']} (Score: {p['priority_score']:.1f})")
            report.append(f"     Категория: {p['category']}")
            report.append(f"     Эффект: {p['expected_effect']}")
            report.append(f"     Риск: {p['risk']}")

        return '\n'.join(report)


def create_enhancement_visualization(analyzer):
    """Создание визуализации стратегий усиления"""
    fig = plt.figure(figsize=(18, 22))

    # Создаем сетку для графиков
    gs = fig.add_gridspec(5, 2, height_ratios=[2, 1.5, 1.5, 2, 2], hspace=0.3, wspace=0.2)

    seq_len = analyzer.length

    # === 1. Карта репрессорных и активаторных сайтов ===
    ax1 = fig.add_subplot(gs[0, :])
    ax1.set_xlim(-50, seq_len + 100)
    ax1.set_ylim(-2, 4)

    # ДНК линия
    ax1.plot([0, seq_len], [0, 0], 'k-', linewidth=4, zorder=1)
    ax1.plot([seq_len, seq_len + 30], [0, 0], 'g-', linewidth=4, zorder=1)

    # ATG маркер
    ax1.annotate('ATG', xy=(seq_len + 15, 0.3), fontsize=12, ha='center',
                fontweight='bold', color='green')

    # Репрессорные сайты (красные, внизу)
    for rep in analyzer.repressor_sites:
        x = rep['start']
        width = rep['end'] - rep['start']
        rect = FancyBboxPatch((x, -1.5), max(width, 8), 0.8,
                              boxstyle="round,pad=0.02", facecolor='#FF6B6B',
                              edgecolor='darkred', linewidth=1.5, alpha=0.8)
        ax1.add_patch(rect)

        # Стрелка "удалить"
        ax1.annotate('✕', xy=(x + width/2, -1.1), fontsize=10,
                    ha='center', va='center', color='darkred', fontweight='bold')

    # Активаторные сайты (зеленые, вверху)
    y_offset = 1
    for act in analyzer.activator_sites:
        x = act['start']
        width = act['end'] - act['start']
        rect = FancyBboxPatch((x, y_offset), max(width, 8), 0.8,
                              boxstyle="round,pad=0.02", facecolor='#90EE90',
                              edgecolor='darkgreen', linewidth=1.5, alpha=0.8)
        ax1.add_patch(rect)

        # Стрелка "усилить"
        ax1.annotate('+', xy=(x + width/2, y_offset + 0.4), fontsize=12,
                    ha='center', va='center', color='darkgreen', fontweight='bold')

    # Потенциальные места для добавления элементов
    potential_positions = [
        (100, 'UAS\nинсерция'),
        (250, 'Gcr1\nдупликация'),
        (seq_len - 30, 'TATA\nоптимизация'),
    ]

    for pos, label in potential_positions:
        ax1.annotate('', xy=(pos, 2.5), xytext=(pos, 3.2),
                    arrowprops=dict(arrowstyle='->', color='blue', lw=2))
        ax1.text(pos, 3.4, label, ha='center', va='bottom', fontsize=8, color='blue')

    # Шкала
    for pos in range(0, seq_len, 100):
        ax1.plot([pos, pos], [-0.2, 0.2], 'k-', linewidth=1)
        ax1.text(pos, -0.4, str(pos+1), ha='center', fontsize=8)

    ax1.set_title('Карта мишеней для модификации промотора', fontsize=14, fontweight='bold')
    ax1.set_ylabel('')
    ax1.set_yticks([])
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['left'].set_visible(False)

    # Легенда
    legend_elements = [
        mpatches.Patch(color='#FF6B6B', label='Репрессорные сайты (удалить)'),
        mpatches.Patch(color='#90EE90', label='Активаторные сайты (усилить)'),
        mpatches.Patch(color='blue', alpha=0.3, label='Потенциальные сайты инсерции')
    ]
    ax1.legend(handles=legend_elements, loc='upper right', fontsize=9)

    # === 2. Диаграмма стратегий снижения репрессии ===
    ax2 = fig.add_subplot(gs[1, 0])

    repression_strategies = [
        ('Делеция CreA\nсайтов', 85, '#FF6B6B'),
        ('Мутация\nCreA сайтов', 70, '#FFA07A'),
        ('creA- штамм', 95, '#FF4500'),
        ('Лимитация\nглюкозы', 60, '#FFB6C1'),
    ]

    names = [s[0] for s in repression_strategies]
    values = [s[1] for s in repression_strategies]
    colors = [s[2] for s in repression_strategies]

    bars = ax2.barh(names, values, color=colors, edgecolor='black', height=0.6)
    ax2.set_xlim(0, 100)
    ax2.set_xlabel('Эффективность снятия репрессии (%)')
    ax2.set_title('Стратегии снижения репрессии', fontsize=12, fontweight='bold')

    for bar, val in zip(bars, values):
        ax2.text(val + 2, bar.get_y() + bar.get_height()/2, f'{val}%',
                va='center', fontsize=10)

    # === 3. Диаграмма стратегий усиления ===
    ax3 = fig.add_subplot(gs[1, 1])

    enhancement_strategies = [
        ('Gcr1 мультипликация', 5, '#90EE90'),
        ('UAS инсерция', 25, '#32CD32'),
        ('TATA оптимизация', 3, '#98FB98'),
        ('Kozak оптимизация', 8, '#00FA9A'),
        ('Гибридный\nпромотор', 50, '#228B22'),
    ]

    names = [s[0] for s in enhancement_strategies]
    values = [s[1] for s in enhancement_strategies]
    colors = [s[2] for s in enhancement_strategies]

    bars = ax3.barh(names, values, color=colors, edgecolor='black', height=0.6)
    ax3.set_xlim(0, 60)
    ax3.set_xlabel('Потенциальное увеличение экспрессии (×)')
    ax3.set_title('Стратегии усиления промотора', fontsize=12, fontweight='bold')

    for bar, val in zip(bars, values):
        ax3.text(val + 1, bar.get_y() + bar.get_height()/2, f'{val}×',
                va='center', fontsize=10)

    # === 4. Схема модификаций CreA сайтов ===
    ax4 = fig.add_subplot(gs[2, :])
    ax4.set_xlim(0, 10)
    ax4.set_ylim(0, 4)
    ax4.axis('off')

    ax4.set_title('Стратегии модификации сайтов CreA/Mig1', fontsize=12, fontweight='bold')

    # Исходный сайт
    ax4.text(0.5, 3.5, 'Исходный:', fontsize=11, fontweight='bold')
    ax4.add_patch(FancyBboxPatch((2, 3.2), 2, 0.6, boxstyle="round",
                                  facecolor='#FFE4E1', edgecolor='red'))
    ax4.text(3, 3.5, 'GCGGGG', fontsize=10, ha='center', family='monospace')
    ax4.text(4.5, 3.5, '→ CreA связывается → Репрессия', fontsize=10)

    # Мутация 1
    ax4.text(0.5, 2.5, 'Мутация 1:', fontsize=11, fontweight='bold')
    ax4.add_patch(FancyBboxPatch((2, 2.2), 2, 0.6, boxstyle="round",
                                  facecolor='#E0FFE0', edgecolor='green'))
    ax4.text(3, 2.5, 'GCAGAG', fontsize=10, ha='center', family='monospace')
    ax4.text(4.5, 2.5, '→ CreA не связывается → Нет репрессии', fontsize=10, color='green')

    # Мутация 2
    ax4.text(0.5, 1.5, 'Мутация 2:', fontsize=11, fontweight='bold')
    ax4.add_patch(FancyBboxPatch((2, 1.2), 2, 0.6, boxstyle="round",
                                  facecolor='#E0FFE0', edgecolor='green'))
    ax4.text(3, 1.5, 'GAGGGG', fontsize=10, ha='center', family='monospace')
    ax4.text(4.5, 1.5, '→ Минимальные изменения GC', fontsize=10, color='green')

    # Делеция
    ax4.text(0.5, 0.5, 'Делеция:', fontsize=11, fontweight='bold')
    ax4.add_patch(FancyBboxPatch((2, 0.2), 2, 0.6, boxstyle="round",
                                  facecolor='#E6E6FA', edgecolor='purple'))
    ax4.text(3, 0.5, '------', fontsize=10, ha='center', family='monospace')
    ax4.text(4.5, 0.5, '→ Полное удаление (проверить спейсинг)', fontsize=10, color='purple')

    # === 5. Схема гибридного промотора ===
    ax5 = fig.add_subplot(gs[3, :])
    ax5.set_xlim(0, 100)
    ax5.set_ylim(0, 8)
    ax5.axis('off')

    ax5.set_title('Конструирование усиленного гибридного промотора', fontsize=12, fontweight='bold')

    # Исходный промотор
    ax5.text(5, 7, 'Исходный pGAP:', fontsize=11, fontweight='bold')
    regions_orig = [
        (10, 25, '#FFE4E1', 'Дистальный'),
        (35, 25, '#E6E6FA', 'Проксим.'),
        (60, 25, '#E0FFE0', 'Коровый'),
        (85, 10, '#FFFACD', 'ATG'),
    ]
    for x, w, c, label in regions_orig:
        ax5.add_patch(FancyBboxPatch((x, 6), w, 0.8, boxstyle="round",
                                      facecolor=c, edgecolor='black'))
        ax5.text(x + w/2, 6.4, label, fontsize=8, ha='center')

    # Модифицированный промотор
    ax5.text(5, 4.5, 'Усиленный:', fontsize=11, fontweight='bold', color='green')

    # UAS кассета
    ax5.add_patch(FancyBboxPatch((10, 3.5), 15, 0.8, boxstyle="round",
                                  facecolor='#90EE90', edgecolor='darkgreen', linewidth=2))
    ax5.text(17.5, 3.9, 'UAS×4', fontsize=9, ha='center', fontweight='bold')

    # Gcr1 кассета
    ax5.add_patch(FancyBboxPatch((27, 3.5), 10, 0.8, boxstyle="round",
                                  facecolor='#98FB98', edgecolor='darkgreen', linewidth=2))
    ax5.text(32, 3.9, 'Gcr1×4', fontsize=9, ha='center', fontweight='bold')

    # Коровый без CreA
    ax5.add_patch(FancyBboxPatch((39, 3.5), 25, 0.8, boxstyle="round",
                                  facecolor='#E0FFE0', edgecolor='green', linewidth=2))
    ax5.text(51.5, 3.9, 'Коровый (ΔCreA)', fontsize=9, ha='center')

    # Оптимизированный TATA
    ax5.add_patch(FancyBboxPatch((66, 3.5), 12, 0.8, boxstyle="round",
                                  facecolor='#FFEAA7', edgecolor='orange', linewidth=2))
    ax5.text(72, 3.9, 'TATAAA', fontsize=9, ha='center', fontweight='bold')

    # Kozak
    ax5.add_patch(FancyBboxPatch((80, 3.5), 10, 0.8, boxstyle="round",
                                  facecolor='#DDA0DD', edgecolor='purple', linewidth=2))
    ax5.text(85, 3.9, 'Kozak', fontsize=9, ha='center', fontweight='bold')

    # ATG
    ax5.add_patch(FancyBboxPatch((92, 3.5), 5, 0.8, boxstyle="round",
                                  facecolor='#00FF00', edgecolor='darkgreen', linewidth=2))
    ax5.text(94.5, 3.9, 'ATG', fontsize=9, ha='center', fontweight='bold')

    # Аннотации улучшений
    improvements = [
        (17.5, 2.8, '↑ 5-10×'),
        (32, 2.8, '↑ 3-5×'),
        (51.5, 2.8, '↑ 2-5×'),
        (72, 2.8, '↑ 2-3×'),
        (85, 2.8, '↑ 2-5×'),
    ]
    for x, y, text in improvements:
        ax5.text(x, y, text, fontsize=9, ha='center', color='green', fontweight='bold')

    # Итоговый эффект
    ax5.add_patch(FancyBboxPatch((25, 1), 50, 1.2, boxstyle="round",
                                  facecolor='#90EE90', edgecolor='darkgreen', linewidth=3))
    ax5.text(50, 1.6, 'Потенциальное увеличение: 50-100× относительно исходного',
            fontsize=11, ha='center', fontweight='bold', color='darkgreen')

    # === 6. Таблица приоритетов ===
    ax6 = fig.add_subplot(gs[4, :])
    ax6.axis('off')

    ax6.set_title('Рекомендованный план модификаций', fontsize=12, fontweight='bold', pad=20)

    priorities = analyzer.calculate_modification_priority()[:8]

    table_data = []
    for i, p in enumerate(priorities, 1):
        table_data.append([
            str(i),
            p['strategy'][:30],
            p['category'].replace('_', ' ').title(),
            p['expected_effect'][:25] + '...' if len(p['expected_effect']) > 25 else p['expected_effect'],
            f"{p['priority_score']:.1f}"
        ])

    table = ax6.table(
        cellText=table_data,
        colLabels=['#', 'Стратегия', 'Категория', 'Эффект', 'Приоритет'],
        loc='center',
        cellLoc='left',
        colWidths=[0.05, 0.3, 0.2, 0.35, 0.1]
    )

    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.8)

    # Окрашивание заголовка
    for j in range(5):
        table[(0, j)].set_facecolor('#4472C4')
        table[(0, j)].set_text_props(color='white', fontweight='bold')

    # Окрашивание по приоритету
    for i, row in enumerate(table_data, 1):
        score = float(row[4])
        if score > 3:
            color = '#C6EFCE'
        elif score > 2:
            color = '#FFEB9C'
        else:
            color = '#FFC7CE'
        for j in range(5):
            table[(i, j)].set_facecolor(color)

    plt.suptitle('АНАЛИЗ СТРАТЕГИЙ УСИЛЕНИЯ ПРОМОТОРА pGAP\nИ СНИЖЕНИЯ РЕПРЕССИИ',
                fontsize=16, fontweight='bold', y=0.98)

    plt.tight_layout(rect=[0, 0, 1, 0.96])

    plt.savefig('promoter_enhancement_strategies.png', dpi=150, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.savefig('promoter_enhancement_strategies.svg', format='svg', bbox_inches='tight')

    print("\n✓ Визуализация сохранена: promoter_enhancement_strategies.png")
    print("✓ Визуализация сохранена: promoter_enhancement_strategies.svg")

    return fig


def create_mutation_map(analyzer):
    """Создание детальной карты рекомендуемых мутаций"""
    fig, ax = plt.subplots(figsize=(16, 10))

    seq = analyzer.sequence
    seq_len = len(seq)

    ax.set_xlim(-20, seq_len + 50)
    ax.set_ylim(-5, 10)

    # Основная линия ДНК
    ax.plot([0, seq_len], [0, 0], 'k-', linewidth=5, zorder=1)

    # Регионы промотора
    regions = [
        (0, seq_len//3, '#FFE4E1', 'Дистальный'),
        (seq_len//3, 2*seq_len//3, '#E6E6FA', 'Проксимальный'),
        (2*seq_len//3, seq_len, '#E0FFE0', 'Коровый'),
    ]

    for start, end, color, name in regions:
        rect = FancyBboxPatch((start, -0.4), end-start, 0.8,
                              boxstyle="round,pad=0.01", facecolor=color,
                              edgecolor='gray', alpha=0.5, zorder=0)
        ax.add_patch(rect)
        ax.text((start+end)/2, 0, name, ha='center', va='center', fontsize=9)

    # Мутации для удаления (репрессоры)
    deletion_y = -2
    for rep in analyzer.repressor_sites:
        x = rep['start']
        width = max(rep['end'] - rep['start'], 10)

        # Красная зона удаления
        rect = FancyBboxPatch((x, deletion_y - 0.4), width, 0.8,
                              boxstyle="round", facecolor='#FF6B6B',
                              edgecolor='darkred', linewidth=2)
        ax.add_patch(rect)

        # Подпись
        ax.text(x + width/2, deletion_y, rep['match'][:10],
               ha='center', va='center', fontsize=7, family='monospace')

        # Линия связи
        ax.plot([x + width/2, x + width/2], [deletion_y + 0.4, -0.4],
               'r--', linewidth=1, alpha=0.5)

        # Метка
        ax.text(x + width/2, deletion_y - 0.8, rep['name'].split()[0],
               ha='center', va='top', fontsize=7, color='darkred')

    # Зоны усиления (активаторы)
    enhance_y = 2
    for act in analyzer.activator_sites:
        x = act['start']
        width = max(act['end'] - act['start'], 10)

        # Зеленая зона усиления
        rect = FancyBboxPatch((x, enhance_y - 0.4), width, 0.8,
                              boxstyle="round", facecolor='#90EE90',
                              edgecolor='darkgreen', linewidth=2)
        ax.add_patch(rect)

        # Подпись
        ax.text(x + width/2, enhance_y, act['match'],
               ha='center', va='center', fontsize=7, family='monospace')

        # Линия связи
        ax.plot([x + width/2, x + width/2], [0.4, enhance_y - 0.4],
               'g--', linewidth=1, alpha=0.5)

        # Метка
        ax.text(x + width/2, enhance_y + 0.8, act['name'].split()[0],
               ha='center', va='bottom', fontsize=7, color='darkgreen')

    # Потенциальные места вставки
    insertion_sites = [
        (50, 'UAS кассета'),
        (seq_len - 80, 'Оптимиз. TATA'),
        (seq_len - 10, 'Kozak'),
    ]

    insert_y = 5
    for pos, label in insertion_sites:
        # Стрелка вставки
        ax.annotate('', xy=(pos, 0.5), xytext=(pos, insert_y - 0.5),
                   arrowprops=dict(arrowstyle='->', color='blue', lw=2))

        # Блок вставки
        rect = FancyBboxPatch((pos - 30, insert_y - 0.4), 60, 0.8,
                              boxstyle="round", facecolor='#87CEEB',
                              edgecolor='darkblue', linewidth=2)
        ax.add_patch(rect)
        ax.text(pos, insert_y, label, ha='center', va='center', fontsize=9)

    # Шкала
    for pos in range(0, seq_len + 1, 100):
        ax.plot([pos, pos], [-0.6, -0.8], 'k-', linewidth=1)
        ax.text(pos, -1.0, f'{pos}', ha='center', fontsize=8)

    # Легенда
    legend_elements = [
        mpatches.Patch(color='#FF6B6B', label='Удалить (репрессоры)'),
        mpatches.Patch(color='#90EE90', label='Усилить (активаторы)'),
        mpatches.Patch(color='#87CEEB', label='Вставить (новые элементы)'),
    ]
    ax.legend(handles=legend_elements, loc='upper left', fontsize=10)

    ax.set_title('Детальная карта рекомендуемых модификаций промотора pGAP',
                fontsize=14, fontweight='bold')
    ax.axis('off')

    plt.tight_layout()
    plt.savefig('promoter_mutation_map.png', dpi=150, bbox_inches='tight')
    plt.savefig('promoter_mutation_map.svg', format='svg', bbox_inches='tight')

    print("✓ Карта мутаций сохранена: promoter_mutation_map.png")

    return fig


def main():
    """Главная функция"""
    print("\n" + "=" * 90)
    print("АНАЛИЗ ВОЗМОЖНОСТЕЙ УСИЛЕНИЯ ПРОМОТОРА pGAP")
    print("И СТРАТЕГИЙ СНИЖЕНИЯ РЕПРЕССИИ")
    print("=" * 90)

    # Создаем анализатор
    analyzer = PromoterEnhancementAnalyzer(PROMOTER_SEQUENCE)

    # Генерируем отчет
    report = analyzer.generate_report()
    print(report)

    # Сохраняем отчет в файл
    with open('promoter_enhancement_report.txt', 'w', encoding='utf-8') as f:
        f.write(report)
    print("\n✓ Отчет сохранен: promoter_enhancement_report.txt")

    # Создаем визуализации
    print("\n" + "=" * 90)
    print("СОЗДАНИЕ ВИЗУАЛИЗАЦИЙ")
    print("=" * 90)

    try:
        fig1 = create_enhancement_visualization(analyzer)
        plt.close(fig1)

        fig2 = create_mutation_map(analyzer)
        plt.close(fig2)

    except Exception as e:
        print(f"Ошибка при создании визуализации: {e}")
        import traceback
        traceback.print_exc()

    # Итоговые рекомендации
    print("\n" + "=" * 90)
    print("КЛЮЧЕВЫЕ РЕКОМЕНДАЦИИ")
    print("=" * 90)
    print("""
╔══════════════════════════════════════════════════════════════════════════════╗
║                    СТРАТЕГИЯ УСИЛЕНИЯ ПРОМОТОРА pGAP                         ║
╠══════════════════════════════════════════════════════════════════════════════╣
║                                                                              ║
║  1. СНИЖЕНИЕ РЕПРЕССИИ (Приоритет: ВЫСОКИЙ)                                  ║
║     ├─ Делеция/мутация сайтов CreA (GCGGGG → GCAGAG)                        ║
║     ├─ Использование creA-дефицитного штамма                                 ║
║     └─ Оптимизация условий культивирования (лимитация глюкозы)               ║
║                                                                              ║
║  2. УСИЛЕНИЕ БАЗАЛЬНОЙ АКТИВНОСТИ (Приоритет: ВЫСОКИЙ)                       ║
║     ├─ Оптимизация TATA-бокса (→ TATAAAA)                                   ║
║     ├─ Оптимизация Kozak (→ ACCATGG)                                        ║
║     └─ Мультипликация Gcr1 сайтов (×4)                                      ║
║                                                                              ║
║  3. ДОБАВЛЕНИЕ UAS ЭЛЕМЕНТОВ (Приоритет: СРЕДНИЙ)                            ║
║     ├─ Вставка UAS из сильных промоторов                                    ║
║     ├─ Тандемные копии энхансеров                                           ║
║     └─ Синтетические активаторные кассеты                                   ║
║                                                                              ║
║  4. КОНСТРУИРОВАНИЕ ГИБРИДНОГО ПРОМОТОРА (Приоритет: СРЕДНИЙ)                ║
║     ├─ UAS(сильный) + Core(pGAP) + оптимизированный TATA                    ║
║     └─ Потенциальное увеличение: 50-100×                                    ║
║                                                                              ║
╚══════════════════════════════════════════════════════════════════════════════╝
    """)

    print("\n" + "=" * 90)
    print("КОНКРЕТНЫЕ МОДИФИКАЦИИ ДЛЯ РЕАЛИЗАЦИИ")
    print("=" * 90)

    crea_sites = [r for r in analyzer.repressor_sites if 'CreA' in r['name']]
    print(f"\nНайдено {len(crea_sites)} сайтов CreA для модификации:")

    for site in crea_sites:
        pos = site['position_from_atg']
        seq = site['match']
        print(f"\n  Позиция {pos:+d}:")
        print(f"    Исходная: 5'-{seq}-3'")

        # Предлагаем мутации
        if 'GCGGGG' in seq:
            print(f"    Мутация:  5'-GCAGAG-3' (минимальные изменения)")
        elif 'GG' in seq:
            mutated = seq.replace('GG', 'AG', 1)
            print(f"    Мутация:  5'-{mutated}-3' (нарушение консенсуса)")

    gcr1_sites = [a for a in analyzer.activator_sites if 'Gcr1' in a['name']]
    print(f"\n\nНайдено {len(gcr1_sites)} сайтов Gcr1 для усиления:")
    print("  Рекомендация: Добавить синтетическую кассету (CTTCC)×4")
    print("  Оптимальная позиция: -300 до -200 от ATG")


if __name__ == "__main__":
    main()
