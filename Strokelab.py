# -*- coding: utf-8 -*-
"""
StrokeLab, этап 3: общая ось буквы и предсказание тонкого начертания.

Разбор по векторным данным шрифта, без растра и без диаграммы Вороного.

Этап 1. Рёбра.
    Контуры глифа раскладываются на рёбра (прямые и дуги); для каждого считаются касательные,
    поворот, длина и нормаль внутрь материала; узлы классифицируются (выпуклый угол, вогнутый
    угол, гладкий узел, перегиб).

Этап 2. Штрихи и их оси.
    У штриха две стороны. Из точки ребра пускается луч по нормали внутрь буквы; пара рёбер
    принимается как две стороны штриха, если нормали смотрят друг на друга и круг с диаметром
    между точками пуст. Штрих строится только между двумя прямыми или двумя дугами. Середины
    отрезков дают ось и профиль толщины. Куски штриха, разорванные перемычкой, склеиваются;
    стыки кусков оси сглаживаются («усы» на одной прямой с узлом, в экстремуме строго по
    горизонтали или вертикали).

Этап 3. Общая ось и предсказание, опираясь на существующий контур.
    * Конец штриха, перед которым контур сразу обрывается, это настоящий торец, соединять его
      не с чем. Если вперёд идёт материал, это стык.
    * Стык принимается только при опоре в контуре: рядом должны быть вогнутые углы (так
      выглядят T, L, X, развилки). Продление идёт только внутри материала; там, где прямая
      выходит из материала (вырез между стойкой и аркой), ось ведётся по середине прохода.
    * Виды соединений: T (конец упирается в ось другого штриха), L (угол двух штрихов),
      Y (три конца в одном узле), X (оси пересекаются), G1 (плавное продолжение: стойка и
      арка, спинка S). Свободный конец достраивается до видимого торца.
    * Предсказание: по осям и измеренной толщине строятся ленты штрихов. При исходной толщине
      они должны ложиться на существующий контур (процент совпадения и «лишнего» показан в
      сведениях), при тонкой толщине дают предсказание тонкого начертания.
    * Два набора настроек: для букв, где все рёбра прямые, и для букв с кривыми. Набор
      выбирается автоматически по классу буквы.
    * Экспорт исходных осевых путей в UFO/Designspace/Glyphs сохраняет открытые линии без заливки.
    * Торцы: лента режется по линии исходного торца, поэтому косые и угловые срезы сохраняются,
      а выступающие за контур углы подрезаются. Точки (i, j, двоеточие) распознаются отдельно и
      в тонком начертании уменьшаются вокруг центра. Стык T по умолчанию прямой: круговой вход
      в стойку применяется только когда конец почти параллелен ей.

Установка (Windows 10, Python 3.9+):
    pip install fonttools PySide6
    pip install skia-pathops      # необязательно: удаление наложений контуров
    pip install glyphsLib ufoLib2  # необязательно: экспорт в Glyphs 3 (UFO/Designspace работают без них)

Запуск:
    python strokelab.py [путь_к_шрифту.otf]
Шрифт можно перетащить мышью в окно.

Управление на холсте:
    колесо мыши       масштаб вокруг курсора
    перетаскивание    сдвиг (левая, средняя или правая кнопка)
    клик              подробности о штрихе (клик по пятну) и ребре (клик у линии)
    двойной клик      вписать глиф в окно
"""
from __future__ import annotations

import math
import os
import sys
import shutil
import tempfile
from types import SimpleNamespace
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from fontTools.pens.basePen import BasePen
from fontTools.ttLib import TTFont

try:  # необязательная зависимость: удаление наложений контуров
    import pathops

    HAVE_PATHOPS = True
except Exception:  # pragma: no cover
    pathops = None
    HAVE_PATHOPS = False

# ============================================================================
#  Часть 1. Геометрия (чистый Python, без Qt)
# ============================================================================

Pt = Tuple[float, float]
Cubic = Tuple[Pt, Pt, Pt, Pt]


def vsub(a: Pt, b: Pt) -> Pt:
    return (a[0] - b[0], a[1] - b[1])


def vdot(a: Pt, b: Pt) -> float:
    return a[0] * b[0] + a[1] * b[1]


def vcross(a: Pt, b: Pt) -> float:
    return a[0] * b[1] - a[1] * b[0]


def dist(a: Pt, b: Pt) -> float:
    return math.hypot(a[0] - b[0], a[1] - b[1])


def vunit(v: Pt) -> Optional[Pt]:
    n = math.hypot(v[0], v[1])
    if n < 1e-9:
        return None
    return (v[0] / n, v[1] / n)


def lerp(a: Pt, b: Pt, t: float) -> Pt:
    return (a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t)


def angle_deg(u: Pt, v: Pt) -> float:
    """Знаковый угол от u к v в градусах, против часовой стрелки положительный
    (система координат шрифта: ось Y вверх)."""
    return math.degrees(math.atan2(vcross(u, v), vdot(u, v)))


def line_cubic(p0: Pt, p3: Pt) -> Cubic:
    """Отрезок, записанный как кубика (контрольные точки на хорде)."""
    return (p0, lerp(p0, p3, 1 / 3), lerp(p0, p3, 2 / 3), p3)


def bez_pt(p: Cubic, t: float) -> Pt:
    mt = 1.0 - t
    a, b, c, d = mt * mt * mt, 3 * mt * mt * t, 3 * mt * t * t, t * t * t
    return (
        a * p[0][0] + b * p[1][0] + c * p[2][0] + d * p[3][0],
        a * p[0][1] + b * p[1][1] + c * p[2][1] + d * p[3][1],
    )


def bez_deriv(p: Cubic, t: float) -> Pt:
    mt = 1.0 - t
    a, b, c = 3 * mt * mt, 6 * mt * t, 3 * t * t
    v1, v2, v3 = vsub(p[1], p[0]), vsub(p[2], p[1]), vsub(p[3], p[2])
    return (a * v1[0] + b * v2[0] + c * v3[0], a * v1[1] + b * v2[1] + c * v3[1])


def split_cubic(p: Cubic, t: float) -> Tuple[Cubic, Cubic]:
    p01, p12, p23 = lerp(p[0], p[1], t), lerp(p[1], p[2], t), lerp(p[2], p[3], t)
    p012, p123 = lerp(p01, p12, t), lerp(p12, p23, t)
    m = lerp(p012, p123, t)
    return (p[0], p01, p012, m), (m, p123, p23, p[3])


def sub_cubic(p: Cubic, t0: float, t1: float) -> Cubic:
    """Кусок кубики на отрезке параметра [t0, t1]."""
    if t1 < 1.0:
        p, _ = split_cubic(p, t1)
    if t0 > 0.0:
        _, p = split_cubic(p, t0 / t1)
    return p


def flatten(p: Cubic, n: int = 16) -> List[Pt]:
    return [bez_pt(p, i / n) for i in range(n + 1)]


def flatten_n(p: Cubic) -> int:
    """Число отрезков для ломаной: примерно по одному на 5° поворота ломаной управляющих точек."""
    vs = [vsub(p[1], p[0]), vsub(p[2], p[1]), vsub(p[3], p[2])]
    turn = 0.0
    for u, v in zip(vs, vs[1:]):
        if vunit(u) and vunit(v):
            turn += abs(angle_deg(u, v))
    return max(12, min(64, int(math.ceil(turn / 5.0))))


def quad_roots(a: float, b: float, c: float) -> List[float]:
    scale = abs(a) + abs(b) + abs(c)
    if scale == 0:
        return []
    if abs(a) <= 1e-9 * scale:
        return [] if abs(b) <= 1e-9 * scale else [-c / b]
    d = b * b - 4 * a * c
    if d < 0:
        return []
    s = math.sqrt(d)
    return [(-b - s) / (2 * a), (-b + s) / (2 * a)]


def critical_ts(p: Cubic, extrema: bool, inflections: bool) -> List[float]:
    """Параметры t внутри кубики, где она должна быть разрезана:
    экстремумы по X и Y и точки перегиба."""
    a, b, c = vsub(p[1], p[0]), vsub(p[2], p[1]), vsub(p[3], p[2])
    # B'(t)/3 = A t^2 + B t + C
    A = (a[0] - 2 * b[0] + c[0], a[1] - 2 * b[1] + c[1])
    B = (2 * (b[0] - a[0]), 2 * (b[1] - a[1]))
    C = a
    ts: List[float] = []
    if extrema:
        ts += quad_roots(A[0], B[0], C[0]) + quad_roots(A[1], B[1], C[1])
    if inflections:
        # cross(B'(t), B''(t)) = -cross(A,B) t^2 + 2 cross(C,A) t + cross(C,B)
        ts += quad_roots(-vcross(A, B), 2 * vcross(C, A), vcross(C, B))
    return sorted(t for t in ts if 1e-4 < t < 1 - 1e-4)


def is_flat(p: Cubic, tol: float) -> bool:
    chord = vsub(p[3], p[0])
    L = math.hypot(chord[0], chord[1])
    if L < 1e-9:
        return False
    d1 = abs(vcross(chord, vsub(p[1], p[0]))) / L
    d2 = abs(vcross(chord, vsub(p[2], p[0]))) / L
    return max(d1, d2) <= tol


def end_tangents(p: Cubic) -> Tuple[Pt, Pt]:
    t0 = vunit(vsub(p[1], p[0])) or vunit(vsub(p[2], p[0])) or vunit(vsub(p[3], p[0]))
    t1 = vunit(vsub(p[3], p[2])) or vunit(vsub(p[3], p[1])) or vunit(vsub(p[3], p[0]))
    return t0 or (1.0, 0.0), t1 or (1.0, 0.0)


def dist_to_polyline(pt: Pt, poly: List[Pt]) -> float:
    best = 1e18
    for a, b in zip(poly, poly[1:]):
        abx, aby = b[0] - a[0], b[1] - a[1]
        l2 = abx * abx + aby * aby
        t = 0.0 if l2 == 0 else max(0.0, min(1.0, ((pt[0] - a[0]) * abx + (pt[1] - a[1]) * aby) / l2))
        d = math.hypot(pt[0] - (a[0] + abx * t), pt[1] - (a[1] + aby * t))
        if d < best:
            best = d
    return best


def winding_number(pt: Pt, polys: List[List[Pt]]) -> int:
    x, y = pt
    wn = 0
    for poly in polys:
        n = len(poly)
        for i in range(n):
            x0, y0 = poly[i]
            x1, y1 = poly[(i + 1) % n]
            is_left = (x1 - x0) * (y - y0) - (x - x0) * (y1 - y0)
            if y0 <= y:
                if y1 > y and is_left > 0:
                    wn += 1
            elif y1 <= y and is_left < 0:
                wn -= 1
    return wn


# ============================================================================
#  Часть 2. Чтение глифов из шрифта
# ============================================================================

Segment = Tuple[str, Tuple[Pt, ...]]  # ('L', (p0, p1)) или ('C', (p0, p1, p2, p3))


class ContourPen(BasePen):
    """Собирает контуры глифа как списки отрезков и кубических кривых.
    Составные глифы раскладываются, квадратичные кривые (TrueType) переводятся
    в кубические средствами BasePen."""

    def __init__(self, glyph_set):
        super().__init__(glyph_set)
        self.contours: List[List[Segment]] = []
        self._cur: List[Segment] = []
        self._start: Pt = (0.0, 0.0)

    def _moveTo(self, pt):
        self._cur = []
        self._start = (float(pt[0]), float(pt[1]))

    def _lineTo(self, pt):
        p0 = self._getCurrentPoint()
        p1 = (float(pt[0]), float(pt[1]))
        if dist(p0, p1) > 1e-6:
            self._cur.append(("L", (p0, p1)))

    def _curveToOne(self, pt1, pt2, pt3):
        p = (
            self._getCurrentPoint(),
            (float(pt1[0]), float(pt1[1])),
            (float(pt2[0]), float(pt2[1])),
            (float(pt3[0]), float(pt3[1])),
        )
        if max(dist(p[0], q) for q in p[1:]) > 1e-6:
            self._cur.append(("C", p))

    def _finish(self):
        cur = self._getCurrentPoint()
        if cur is not None and dist(cur, self._start) > 1e-6:
            self._cur.append(("L", (cur, self._start)))
        if self._cur:
            self.contours.append(self._cur)
        self._cur = []

    def _closePath(self):
        self._finish()

    def _endPath(self):
        self._finish()


@dataclass
class RawGlyph:
    name: str
    contours: List[List[Segment]]
    advance: float

    def bbox(self) -> Optional[Tuple[float, float, float, float]]:
        xs, ys = [], []
        for c in self.contours:
            for _, pts in c:
                for p in pts:
                    xs.append(p[0])
                    ys.append(p[1])
        if not xs:
            return None
        return min(xs), min(ys), max(xs), max(ys)


class FontStore:
    """Шрифт и кэш прочитанных глифов."""

    def __init__(self, path: str):
        self.path = path
        self.font = TTFont(path, lazy=True)
        self.glyph_set = self.font.getGlyphSet()
        self.upm = int(self.font["head"].unitsPerEm)
        try:
            self.ascent = float(self.font["hhea"].ascent)
            self.descent = float(self.font["hhea"].descent)
        except Exception:
            self.ascent, self.descent = 0.8 * self.upm, -0.2 * self.upm
        try:
            self.title = self.font["name"].getDebugName(4) or os.path.basename(path)
        except Exception:
            self.title = os.path.basename(path)

        cmap = self.font.getBestCmap() or {}
        self.char_of: Dict[str, str] = {}
        self.cp_of: Dict[str, int] = {}
        ordered: List[str] = []
        for cp in sorted(cmap):
            name = cmap[cp]
            if name not in self.char_of and name in self.glyph_set:
                self.char_of[name] = chr(cp)
                self.cp_of[name] = cp
                ordered.append(name)
        seen = set(ordered)
        for name in self.font.getGlyphOrder():
            if name not in seen:
                ordered.append(name)
        self.display_order = ordered
        self.name_of_char = {c: n for n, c in self.char_of.items()}
        self._cache: Dict[Tuple[str, bool], RawGlyph] = {}

    def raw(self, name: str, remove_overlaps: bool = False) -> RawGlyph:
        key = (name, bool(remove_overlaps and HAVE_PATHOPS))
        r = self._cache.get(key)
        if r is None:
            r = self._load(name, key[1])
            self._cache[key] = r
        return r

    def _load(self, name: str, remove_overlaps: bool) -> RawGlyph:
        g = self.glyph_set[name]
        pen = ContourPen(self.glyph_set)
        done = False
        if remove_overlaps:
            try:
                path = pathops.Path()
                g.draw(path.getPen(glyphSet=self.glyph_set))
                path.simplify(fix_winding=True, keep_starting_points=False)
                path.draw(pen)
                done = True
            except Exception:
                pen = ContourPen(self.glyph_set)
        if not done:
            g.draw(pen)
        return RawGlyph(name, pen.contours, float(getattr(g, "width", 0) or 0))

    def describe(self, name: str) -> str:
        ch = self.char_of.get(name)
        if ch:
            return f"{name}  (U+{ord(ch):04X} «{ch}»)"
        return name


# ============================================================================
#  Часть 3. Рёбра, узлы, нормали
# ============================================================================


@dataclass
class Params:
    smooth_deg: float = 45.0  # узел с изломом не больше этого угла считается гладким
    flat_tol: float = 4.0  # кривая, отклоняющаяся от хорды не больше, считается прямой (ед. шрифта)
    split_extrema: bool = False  # резать кривые в экстремумах
    split_inflections: bool = False  # резать кривые в точках перегиба
    remove_overlaps: bool = True  # удалять наложения контуров (нужен skia-pathops)
    max_dev: float = 22.0  # макс. отклонение нормалей двух сторон штриха, градусы
    disc_slack: float = 1.0  # строгость проверки пустого круга (0.5 мягко .. 1.0 строго)
    merge_deg: float = 45.0  # допуск угла при склейке осей-дуг, градусы
    line_merge_deg: float = 5.0  # допуск угла при склейке осей-прямых, градусы (меньше, иначе ветвь Y слипается со стойкой)
    gap_factor: float = 0.5  # макс. разрыв при склейке, в толщинах штриха
    blob_ratio: float = 0.5  # штрих короче этого числа толщин считается пятном (0: не применять)
    allow_mixed: bool = True  # разрешить пары «прямая + дуга» (по умолчанию такие пары считаются ошибочными)
    dot_ratio: float = 1.5  # отдельная деталь не крупнее этого числа толщин: точка (i, j, двоеточие); 0: не искать
    smooth_joints: bool = True  # сглаживать стыки кусков оси (усы на одной прямой с узлом)
    snap_deg: float = 25.0  # касательная в стыке ближе этого угла к горизонтали или вертикали: точно по ней
    join_axes: bool = True  # соединять оси штрихов в общую ось
    ext_factor: float = 1.2  # макс. продление оси к стыку, в толщинах штриха
    min_support: int = 2  # сколько вогнутых углов контура должно быть рядом со стыком
    term_factor: float = 0.5  # конец штриха, перед которым торец ближе этой доли толщины, считается торцом
    cap_factor: float = 5.0  # свободный конец достраивается до видимого торца, если он ближе (в толщинах); 0: не достраивать
    thin_ratio: float = 0.4  # толщина предсказываемого тонкого начертания: доля от исходной
    arm_continue: bool = True  # чаша b, d, p, q: от вершины чаши ось идёт одной дугой окружности и садится на ось стойки под углом arm_angle
    arm_angle: float = 38.0  # угол входа этой дуги в ось стойки, градусы (меньше: положе и ниже в стойке, больше: круче и выше)
    medial_joints: bool = True  # стыки дуг: ось ведётся по настоящей середине контура (круг, касающийся двух сторон), а не по шаблонной кривой


@dataclass(eq=False)
class Edge:
    contour: int
    index: int  # номер ребра внутри контура
    gid: int  # сквозной номер ребра в глифе
    kind: str  # 'line' | 'curve'
    pts: Cubic  # для прямой контрольные точки лежат на хорде
    length: float
    t_in: Pt  # единичная касательная в начале
    t_out: Pt  # единичная касательная в конце
    turn: float  # полный поворот касательной, градусы, + против часовой
    poly: List[Pt]
    bbox: Tuple[float, float, float, float] = (0.0, 0.0, 0.0, 0.0)
    cover: List[Tuple[float, float]] = field(default_factory=list, repr=False)
    fill_left: bool = True  # материал слева от направления обхода?
    bulge: str = "straight"  # 'straight' | 'convex' | 'concave'
    v_start: Optional["Vertex"] = field(default=None, repr=False)
    v_end: Optional["Vertex"] = field(default=None, repr=False)

    def point(self, t: float) -> Pt:
        return bez_pt(self.pts, t)

    def uncovered(self, min_len: float = 2.0) -> List[Tuple[float, float]]:
        """Участки ребра (в долях t), не вошедшие ни в один штрих."""
        res, cur = [], 0.0
        for a, b in sorted(self.cover):
            if a > cur:
                res.append((cur, a))
            cur = max(cur, b)
        if cur < 1.0:
            res.append((cur, 1.0))
        return [(a, b) for a, b in res if (b - a) * self.length >= min_len]

    def tangent(self, t: float) -> Pt:
        if self.kind == "line":
            return self.t_in
        u = vunit(bez_deriv(self.pts, t))
        return u or (self.t_in if t < 0.5 else self.t_out)

    def normal(self, t: float) -> Pt:
        """Единичная нормаль, направленная внутрь материала."""
        tx, ty = self.tangent(t)
        return (-ty, tx) if self.fill_left else (ty, -tx)


@dataclass(eq=False)
class Vertex:
    pos: Pt
    angle: float  # излом касательной, градусы
    kind: str  # 'convex' | 'concave' | 'smooth' | 'inflection'
    contour: int
    edge_in: Edge = field(repr=False)
    edge_out: Edge = field(repr=False)


@dataclass
class Analysis:
    contours: List[List[Edge]]
    vertices: List[Vertex]
    raw_contours: List[List[Segment]]
    advance: float
    bbox: Tuple[float, float, float, float]
    edges: List[Edge] = field(default_factory=list)
    polys: List[List[Pt]] = field(default_factory=list)
    strokes: List["Stroke"] = field(default_factory=list)
    stroke_error: str = ""
    glyph_class: str = "curved"  # 'straight' (только прямые рёбра) | 'curved'
    joints: List["Joint"] = field(default_factory=list)
    end_marks: List[Tuple[Pt, str]] = field(default_factory=list)
    pred_orig: List[List[Pt]] = field(default_factory=list)  # ленты штрихов при исходной толщине
    pred_thin: List[List[Pt]] = field(default_factory=list)  # ленты при тонкой толщине
    scores: Dict[str, float] = field(default_factory=dict)
    skeleton_error: str = ""
    cn_strokes: List["Stroke"] = field(default_factory=list, repr=False)  # штрихи с осями (без точек и пятен)
    cn_corners: list = field(default_factory=list, repr=False)  # угловые стыки для среза митры
    mean_w: float = 0.0

    def counts(self) -> Dict[str, int]:
        c = {"contours": len(self.contours), "edges": len(self.edges)}
        c["lines"] = sum(1 for e in self.edges if e.kind == "line")
        c["curves"] = c["edges"] - c["lines"]
        for k in ("convex", "concave", "smooth", "inflection"):
            c[k] = sum(1 for v in self.vertices if v.kind == k)
        for k in ("line", "arc", "blob", "dot"):
            c["s_" + k] = sum(1 for s in self.strokes if s.kind == k)
        for k in ("T", "L", "Y", "X", "G1"):
            c["j_" + k] = sum(1 for j in self.joints if j.kind == k)
        c["open_ends"] = sum(1 for _, m in self.end_marks if m == "open")
        return c


def _merge_lines(items: List[Tuple[str, Cubic]], tol_deg: float = 0.5) -> List[Tuple[str, Cubic]]:
    """Склеивает соседние почти коллинеарные прямые (лишние узлы на прямой)."""
    sin_tol = math.sin(math.radians(tol_deg))

    def mergeable(a, b):
        if a[0] != "line" or b[0] != "line":
            return False
        u = vunit(vsub(a[1][3], a[1][0]))
        v = vunit(vsub(b[1][3], b[1][0]))
        return bool(u and v and vdot(u, v) > 0 and abs(vcross(u, v)) < sin_tol)

    res: List[Tuple[str, Cubic]] = []
    for it in items:
        if res and mergeable(res[-1], it):
            res[-1] = ("line", line_cubic(res[-1][1][0], it[1][3]))
        else:
            res.append(it)
    if len(res) > 1 and mergeable(res[-1], res[0]):
        res[0] = ("line", line_cubic(res[-1][1][0], res[0][1][3]))
        res.pop()
    return res


def _contour_items(segs: List[Segment], prm: Params, min_len: float) -> List[Tuple[str, Cubic]]:
    items: List[Tuple[str, Cubic]] = []
    for kind, pts in segs:
        if kind == "L":
            items.append(("line", line_cubic(pts[0], pts[1])))
            continue
        p: Cubic = pts  # type: ignore[assignment]
        if is_flat(p, prm.flat_tol):
            items.append(("line", line_cubic(p[0], p[3])))
            continue
        bounds = [0.0]
        for t in critical_ts(p, prm.split_extrema, prm.split_inflections):
            q = bez_pt(p, t)
            if dist(q, bez_pt(p, bounds[-1])) < min_len or dist(q, p[3]) < min_len:
                continue
            bounds.append(t)
        bounds.append(1.0)
        for a, b in zip(bounds, bounds[1:]):
            piece = sub_cubic(p, a, b)
            if max(dist(piece[0], q) for q in piece[1:]) < 1e-6:
                continue
            items.append(("curve", piece))
    return _merge_lines(items)


def _make_edge(kind: str, pts: Cubic, ci: int, idx: int) -> Edge:
    t_in, t_out = end_tangents(pts)
    if kind == "line":
        poly = [pts[0], pts[3]]
        length = dist(pts[0], pts[3])
        turn = 0.0
    else:
        poly = flatten(pts, flatten_n(pts))
        length = sum(dist(a, b) for a, b in zip(poly, poly[1:]))
        turn = angle_deg(t_in, t_out)
    xs, ys = [p[0] for p in poly], [p[1] for p in poly]
    bbox = (min(xs) - 0.5, min(ys) - 0.5, max(xs) + 0.5, max(ys) + 0.5)
    return Edge(ci, idx, -1, kind, pts, length, t_in, t_out, turn, poly, bbox)


def _fill_is_left(edges: List[Edge], polys: List[List[Pt]], eps: float) -> bool:
    """Определяет, с какой стороны от направления обхода лежит материал.
    Проверка по правилу ненулевого числа оборотов на трёх самых длинных рёбрах."""
    score = 0
    for e in sorted(edges, key=lambda e: -e.length)[:3]:
        p = e.point(0.5)
        tx, ty = e.tangent(0.5)
        nl = (-ty, tx)
        left = winding_number((p[0] + nl[0] * eps, p[1] + nl[1] * eps), polys) != 0
        right = winding_number((p[0] - nl[0] * eps, p[1] - nl[1] * eps), polys) != 0
        if left and not right:
            score += 1
        elif right and not left:
            score -= 1
    return score >= 0


def analyze(raw: RawGlyph, upm: int, profiles) -> Analysis:
    """profiles: один Params или словарь {'straight': Params, 'curved': Params}.
    Рёбра строятся по настройкам набора «curved» (параметры рёбер общие для обоих наборов),
    затем по классу буквы выбирается набор для поиска штрихов и соединения осей."""
    if not isinstance(profiles, dict):
        profiles = {"straight": profiles, "curved": profiles}
    prm = profiles["curved"]
    min_len = max(0.3, upm * 0.0005)
    contours: List[List[Edge]] = []
    for segs in raw.contours:
        items = _contour_items(segs, prm, min_len)
        edges = [_make_edge(k, p, len(contours), i) for i, (k, p) in enumerate(items)]
        if edges:
            contours.append(edges)

    polys: List[List[Pt]] = []
    for edges in contours:
        poly: List[Pt] = []
        for e in edges:
            poly.extend(e.poly[:-1])
        polys.append(poly)

    eps = max(0.5, upm * 0.001)
    vertices: List[Vertex] = []
    gid = 0
    for ci, edges in enumerate(contours):
        fill_left = _fill_is_left(edges, polys, eps)
        side = 1.0 if fill_left else -1.0
        for e in edges:
            e.fill_left = fill_left
            e.gid = gid
            gid += 1
            if e.kind == "curve":
                e.bulge = "convex" if e.turn * side > 0 else "concave"
        n = len(edges)
        for i in range(n):
            e_in, e_out = edges[i - 1], edges[i]
            ang = angle_deg(e_in.t_out, e_out.t_in)
            if abs(ang) <= prm.smooth_deg:
                if (
                    e_in.kind == "curve"
                    and e_out.kind == "curve"
                    and e_in.turn * e_out.turn < 0
                    and min(abs(e_in.turn), abs(e_out.turn)) > 1.0
                ):
                    kind = "inflection"
                else:
                    kind = "smooth"
            else:
                kind = "convex" if ang * side > 0 else "concave"
            v = Vertex(e_out.pts[0], ang, kind, ci, e_in, e_out)
            e_in.v_end = v
            e_out.v_start = v
            vertices.append(v)

    all_edges = [e for c in contours for e in c]
    bb = raw.bbox() or (0.0, 0.0, raw.advance, 0.0)
    res = Analysis(contours, vertices, raw.contours, raw.advance, bb, all_edges, polys)
    res.glyph_class = "straight" if all(e.kind == "line" for e in all_edges) else "curved"
    sp = profiles[res.glyph_class]
    try:
        find_strokes(res, upm, sp)
    except Exception as ex:  # разбор на рёбра остаётся доступным при сбое поиска штрихов
        res.stroke_error = f"{type(ex).__name__}: {ex}"
        res.strokes = []
    try:
        connect_axes(res, upm, sp)
    except Exception as ex:
        res.skeleton_error = f"{type(ex).__name__}: {ex}"
        res.joints, res.end_marks, res.pred_orig, res.pred_thin, res.scores = [], [], [], [], {}
    return res


# ============================================================================
#  Часть 3б. Штрихи: пары противоположных рёбер, оси, склейка
# ============================================================================
#
#  Идея. У штриха две стороны. Из точки A(t) ребра A пускаем луч вдоль нормали
#  внутрь материала. Первое пересечение с контуром даёт точку B(s) на другом
#  ребре. Пара (A, B) годится как стороны штриха, если:
#    1. нормали сторон смотрят навстречу друг другу (угол не больше max_dev);
#    2. круг с диаметром A(t)-B(s) пуст: другие рёбра в него не заходят.
#       Это отсекает ложные пары через перемычки (стороны H через перекладину);
#    3. таких лучей подряд достаточно много (сторона не просто задета лучом).
#  Середины отрезков A(t)-B(s) образуют ось, длины дают профиль толщины.
#  Куски штриха, разорванные перемычкой (стойка H, перекладина T, диагонали x),
#  склеиваются, если их оси продолжают друг друга.

Sample = Tuple[float, float, Pt, Pt, float]  # (t на A, t на B, точка на A, точка на B, толщина)


@dataclass(eq=False)
class Piece:
    """Кусок штриха: участок стороны A и участок стороны B напротив него."""

    a: Edge
    b: Edge
    samples: List[Sample]
    kind: str  # 'line' | 'arc' | 'blob'
    axis: Cubic
    width: float  # средняя толщина
    length: float  # длина оси
    t_start: Pt
    t_end: Pt
    ta: Tuple[float, float]
    tb: Tuple[float, float]
    axis_parts: List[Cubic] = field(default_factory=list, repr=False)


@dataclass(eq=False)
class Stroke:
    sid: int
    kind: str  # 'line' (прямой) | 'arc' (дуга) | 'blob' (пятно без оси) | 'dot' (точка)
    pieces: List[Piece]
    axis: List[Cubic]  # пусто для пятна
    patches: List[List[Pt]]  # многоугольники для показа
    w_min: float
    w_max: float
    w_mean: float
    length: float
    closed: bool = False
    hue: float = 0.0
    contour: int = -1  # для точки: номер её контура
    caps: Dict[int, Tuple[Pt, Pt]] = field(default_factory=dict, repr=False)  # 0/1 -> (точка, направление торца)
    tags: List[str] = field(default_factory=list)  # по сегменту оси: 'fit' | 'bridge' | 'ext' | 'cap'
    prof: List[Tuple[Pt, float]] = field(default_factory=list, repr=False)  # (середина A-B, толщина)

    def edge_ids(self) -> List[int]:
        ids = set()
        for p in self.pieces:
            ids.add(p.a.gid)
            ids.add(p.b.gid)
        return sorted(ids)


def fit_axis(m: List[Pt], flat_tol: float) -> Tuple[str, Cubic]:
    """Аппроксимирует измеренные середины отрезков прямой или кубикой Безье.

    В отличие от прежней версии, короткий набор из трёх точек не превращается
    автоматически в прямую. Для кривых используется параметризация по длине
    измеренной ломаной и аппроксимация с учётом касательных на концах.
    """
    if not m:
        p = (0.0, 0.0)
        return "line", line_cubic(p, p)

    # Удаляем подряд идущие практически совпадающие измерения: они ухудшают
    # параметризацию и могут сделать систему аппроксимации вырожденной.
    pts: List[Pt] = [m[0]]
    for q in m[1:]:
        if dist(pts[-1], q) > 1e-7:
            pts.append(q)
    m = pts
    if len(m) == 1:
        return "line", line_cubic(m[0], m[0])

    p0, p3 = m[0], m[-1]
    chord_vec = vsub(p3, p0)
    chord = math.hypot(chord_vec[0], chord_vec[1])
    if chord < 1e-6:
        # Замкнувшийся короткий фрагмент не имеет надёжного направления.
        return "arc", (p0, m[len(m) // 3], m[(2 * len(m)) // 3], p3)
    if len(m) == 2:
        return "line", line_cubic(p0, p3)

    dev = max(abs(vcross(chord_vec, vsub(q, p0))) / chord for q in m)
    if dev <= max(flat_tol, 0.01 * chord):
        return "line", line_cubic(p0, p3)

    # Параметры точек по длине ломаной, а не по равномерному индексу.
    d = [0.0]
    for i in range(1, len(m)):
        d.append(d[-1] + dist(m[i - 1], m[i]))
    total = d[-1]
    if total < 1e-8:
        return "line", line_cubic(p0, p3)
    ts = [x / total for x in d]

    def max_fit_error(c: Cubic) -> float:
        return max(dist(bez_pt(c, t), q) for t, q in zip(ts, m))

    # При трёх точках строим квадратичную кривую, преобразованную в кубическую.
    # Она точно проходит через среднее измерение и не получает неустойчивые ручки.
    if len(m) == 3:
        t = ts[1]
        den = 2.0 * (1.0 - t) * t
        if den > 1e-9:
            mt = 1.0 - t
            qmid = ((m[1][0] - mt * mt * p0[0] - t * t * p3[0]) / den,
                    (m[1][1] - mt * mt * p0[1] - t * t * p3[1]) / den)
            q1 = (p0[0] + (2.0 / 3.0) * (qmid[0] - p0[0]),
                  p0[1] + (2.0 / 3.0) * (qmid[1] - p0[1]))
            q2 = (p3[0] + (2.0 / 3.0) * (qmid[0] - p3[0]),
                  p3[1] + (2.0 / 3.0) * (qmid[1] - p3[1]))
            return "arc", (p0, q1, q2, p3)

    # Касательные берутся по крайним измеренным участкам. Это сохраняет
    # направление оси в местах, где фрагмент заканчивается у перемычки.
    t0 = vunit(vsub(m[1], m[0])) or vunit(chord_vec) or (1.0, 0.0)
    t1 = vunit(vsub(m[-1], m[-2])) or vunit(chord_vec) or (1.0, 0.0)

    # Сначала строим кубику с управляющими точками на касательных:
    # P1 = P0 + alpha*T0, P2 = P3 - beta*T1.
    # alpha и beta подбираются совместно методом наименьших квадратов.
    a11 = a12 = a22 = r1 = r2 = 0.0
    tangent_dot = vdot(t0, t1)
    for t, q in zip(ts[1:-1], m[1:-1]):
        mt = 1.0 - t
        b0 = mt * mt * mt
        b1 = 3.0 * mt * mt * t
        b2 = 3.0 * mt * t * t
        b3 = t * t * t
        base = ((b0 + b1) * p0[0] + (b2 + b3) * p3[0],
                (b0 + b1) * p0[1] + (b2 + b3) * p3[1])
        residual = vsub(q, base)
        a11 += b1 * b1
        a12 -= b1 * b2 * tangent_dot
        a22 += b2 * b2
        r1 += b1 * vdot(t0, residual)
        r2 -= b2 * vdot(t1, residual)

    det = a11 * a22 - a12 * a12
    constrained: Optional[Cubic] = None
    if abs(det) > 1e-12:
        alpha = (r1 * a22 - r2 * a12) / det
        beta = (a11 * r2 - a12 * r1) / det
        # Не допускаем отрицательные или чрезмерно длинные ручки: они часто
        # означают, что несколько точек описывают сложную форму, а не одну дугу.
        handle_limit = 2.0 * total
        alpha = max(0.0, min(handle_limit, alpha))
        beta = max(0.0, min(handle_limit, beta))
        p1 = (p0[0] + t0[0] * alpha, p0[1] + t0[1] * alpha)
        p2 = (p3[0] - t1[0] * beta, p3[1] - t1[1] * beta)
        constrained = (p0, p1, p2, p3)

    # Обычная свободная кубическая аппроксимация остаётся запасным вариантом
    # для кривых, чьи крайние касательные оценены нестабильно.
    free: Optional[Cubic] = None
    if len(m) >= 4:
        f11 = f12 = f22 = 0.0
        fr1 = [0.0, 0.0]
        fr2 = [0.0, 0.0]
        for t, q in zip(ts, m):
            mt = 1.0 - t
            b0, b1, b2, b3 = mt ** 3, 3 * mt * mt * t, 3 * mt * t * t, t ** 3
            f11 += b1 * b1
            f12 += b1 * b2
            f22 += b2 * b2
            for k in (0, 1):
                residual = q[k] - b0 * p0[k] - b3 * p3[k]
                fr1[k] += b1 * residual
                fr2[k] += b2 * residual
        fdet = f11 * f22 - f12 * f12
        if abs(fdet) > 1e-12:
            fp1 = tuple((f22 * fr1[k] - f12 * fr2[k]) / fdet for k in (0, 1))
            fp2 = tuple((f11 * fr2[k] - f12 * fr1[k]) / fdet for k in (0, 1))
            free = (p0, fp1, fp2, p3)  # type: ignore[assignment]

    candidates = [c for c in (constrained, free) if c is not None]
    if not candidates:
        return "arc", (p0, lerp(p0, p3, 1 / 3), lerp(p0, p3, 2 / 3), p3)

    # Предпочитаем касательно-согласованную кривую, пока её ошибка не намного
    # больше ошибки свободной аппроксимации. Это уменьшает изломы между частями.
    errors = {id(c): max_fit_error(c) for c in candidates}
    best = min(candidates, key=lambda c: errors[id(c)])
    if constrained is not None and best is not constrained:
        constrained_error = errors[id(constrained)]
        best_error = errors[id(best)]
        if constrained_error <= best_error * 1.15 + max(0.25, flat_tol * 0.25):
            best = constrained
    return "arc", best


def fit_axis_segments(m: List[Pt], flat_tol: float, max_error: float,
                      max_depth: int = 4) -> List[Cubic]:
    """Адаптивно аппроксимирует срединную линию одной или несколькими кубиками.

    Сложный изгиб не заставляется помещаться в одну кубику: если максимальное
    отклонение от измеренных середин слишком велико, набор делится в наиболее
    проблемной внутренней точке. Соседние кубики при этом имеют общий узел.
    """
    if len(m) <= 4:
        return [fit_axis(m, flat_tol)[1]]

    kind, curve = fit_axis(m, flat_tol)
    if kind == "line":
        return [curve]

    distances = [0.0]
    for a, b in zip(m, m[1:]):
        distances.append(distances[-1] + dist(a, b))
    total = distances[-1]
    if total <= 1e-8:
        return [curve]
    params = [d / total for d in distances]
    errors = [dist(bez_pt(curve, t), q) for t, q in zip(params, m)]
    worst = max(range(1, len(m) - 1), key=lambda i: errors[i])
    if errors[worst] <= max_error or max_depth <= 0:
        return [curve]

    # Не допускаем вырожденные фрагменты из двух точек: обе части должны
    # содержать хотя бы три измерения, если исходных данных достаточно.
    split = max(2, min(len(m) - 3, worst))
    if split < 2 or split > len(m) - 3:
        return [curve]
    left = fit_axis_segments(m[:split + 1], flat_tol, max_error, max_depth - 1)
    right = fit_axis_segments(m[split:], flat_tol, max_error, max_depth - 1)
    return left + right


def close_gaps(chain: List[Cubic], closed: bool, tol: float) -> List[Cubic]:
    """Соседние кубики оси, концы которых разошлись не больше чем на tol, сводятся в общий
    узел (середина щели); усы сдвигаются вместе с концами."""
    n = len(chain)
    if n < 2:
        return list(chain)
    segs = [list(c) for c in chain]
    for j in (range(n) if closed else range(1, n)):
        a, b = segs[(j - 1) % n], segs[j % n]
        if 1e-9 < dist(a[3], b[0]) <= tol:
            node = lerp(a[3], b[0], 0.5)
            da = vsub(node, a[3])
            db = vsub(node, b[0])
            a[2] = (a[2][0] + da[0], a[2][1] + da[1])
            a[3] = node
            b[1] = (b[1][0] + db[0], b[1][1] + db[1])
            b[0] = node
    return [tuple(c) for c in segs]  # type: ignore[misc]


def smooth_chain(chain: List[Cubic], closed: bool, snap_deg: float) -> List[Cubic]:
    """Сглаживает стыки соседних кубик оси: «усы» по обе стороны узла ложатся на одну
    прямую с узлом (длины усов сохраняются). Если эта касательная ближе snap_deg к
    горизонтали или вертикали, то есть узел стоит в экстремуме, она берётся точно по ним."""
    n = len(chain)
    if n < 2:
        return list(chain)
    segs = [list(c) for c in chain]
    for j in (range(n) if closed else range(1, n)):
        a, b = segs[(j - 1) % n], segs[j % n]
        ta = end_tangents(tuple(a))[1]  # type: ignore[arg-type]
        tb = end_tangents(tuple(b))[0]  # type: ignore[arg-type]
        t = vunit((ta[0] + tb[0], ta[1] + tb[1])) or ta
        if snap_deg > 0:
            ang = math.degrees(math.atan2(t[1], t[0]))
            for target in (0.0, 90.0, 180.0, -90.0):
                if abs((ang - target + 180.0) % 360.0 - 180.0) <= snap_deg:
                    r = math.radians(target)
                    t = (math.cos(r), math.sin(r))
                    break
        la = dist(a[2], a[3]) or dist(a[0], a[3]) / 3.0
        lb = dist(b[0], b[1]) or dist(b[0], b[3]) / 3.0
        a[2] = (a[3][0] - t[0] * la, a[3][1] - t[1] * la)
        b[1] = (b[0][0] + t[0] * lb, b[0][1] + t[1] * lb)
    return [tuple(c) for c in segs]  # type: ignore[misc]


def convex_hull(points: List[Pt]) -> List[Pt]:
    pts = sorted(set(points))
    if len(pts) <= 2:
        return pts

    def half(seq):
        h: List[Pt] = []
        for p in seq:
            while len(h) >= 2 and vcross(vsub(h[-1], h[-2]), vsub(p, h[-2])) <= 0:
                h.pop()
            h.append(p)
        return h

    lo, up = half(pts), half(reversed(pts))
    return lo[:-1] + up[:-1]


def _reverse_cubic(c: Cubic) -> Cubic:
    return (c[3], c[2], c[1], c[0])


def _cubic_len(c: Cubic) -> float:
    f = flatten(c, 12)
    return sum(dist(a, b) for a, b in zip(f, f[1:]))


def _ray_may_hit(o: Pt, d: Pt, bb) -> bool:
    tmin, tmax = 0.0, 1e18
    for oo, dd, lo, hi in ((o[0], d[0], bb[0], bb[2]), (o[1], d[1], bb[1], bb[3])):
        if abs(dd) < 1e-12:
            if oo < lo or oo > hi:
                return False
        else:
            ta, tb = (lo - oo) / dd, (hi - oo) / dd
            if ta > tb:
                ta, tb = tb, ta
            if ta > tmin:
                tmin = ta
            if tb < tmax:
                tmax = tb
            if tmin > tmax:
                return False
    return True


def _ray_edge(e: Edge, o: Pt, d: Pt, eps: float) -> Optional[Tuple[float, float]]:
    """Ближайшее пересечение луча o + s*d (s > eps) с ребром: (s, t)."""
    poly = e.poly
    n = len(poly) - 1
    best = None
    for i in range(n):
        p0, p1 = poly[i], poly[i + 1]
        ex, ey = p1[0] - p0[0], p1[1] - p0[1]
        den = d[0] * ey - d[1] * ex
        if abs(den) < 1e-12:
            continue
        wx, wy = p0[0] - o[0], p0[1] - o[1]
        s = (wx * ey - wy * ex) / den
        u = (wx * d[1] - wy * d[0]) / den
        if s > eps and -1e-7 <= u <= 1 + 1e-7 and (best is None or s < best[0]):
            best = (s, (i + u) / n)
    if best is None:
        return None
    s, t = best
    if e.kind == "curve":  # уточняем по настоящей кубике (Ньютон)
        for _ in range(4):
            q = bez_pt(e.pts, t)
            dq = bez_deriv(e.pts, t)
            g = (q[0] - o[0]) * d[1] - (q[1] - o[1]) * d[0]
            gp = dq[0] * d[1] - dq[1] * d[0]
            if abs(gp) < 1e-12:
                break
            t = min(1.0, max(0.0, t - g / gp))
        q = bez_pt(e.pts, t)
        s = (q[0] - o[0]) * d[0] + (q[1] - o[1]) * d[1]
        if s <= eps:
            return None
    return s, t


class _StrokeFinder:
    def __init__(self, a: "Analysis", upm: int, prm: Params):
        self.a, self.upm, self.prm = a, upm, prm
        self.edges = a.edges
        self.cos_dev = math.cos(math.radians(prm.max_dev))
        self.eps = max(0.05, upm * 0.00005)
        self._excl: Dict[Tuple[int, int], set] = {}
        self.w_est = self.estimate_width()

    def estimate_width(self) -> float:
        """Типичная толщина штриха глифа: нижняя квартиль первых попаданий лучей по нормалям.
        От неё зависит шаг выборки: короткие детали строчных (рука f, перекладина t) должны
        получать столько же лучей, сколько длинные детали прописных."""
        ws: List[float] = []
        for A in self.edges:
            for t in (0.25, 0.5, 0.75):
                o, d = A.point(t), A.normal(t)
                best = None
                for e in self.edges:
                    if e is A or not _ray_may_hit(o, d, e.bbox):
                        continue
                    h = _ray_edge(e, o, d, self.eps)
                    if h is not None and (best is None or h[0] < best[0]):
                        best = (h[0], e, h[1])
                if best is not None:
                    nB = best[1].normal(best[2])
                    if -(d[0] * nB[0] + d[1] * nB[1]) >= self.cos_dev:
                        ws.append(best[0])
        if not ws:
            return self.upm * 0.1
        ws.sort()
        return ws[len(ws) // 4]

    # --- лучи ------------------------------------------------------------------
    @staticmethod
    def _neighbors(e: Edge) -> set:
        res = {e.gid}
        if e.v_start:
            res.add(e.v_start.edge_in.gid)
        if e.v_end:
            res.add(e.v_end.edge_out.gid)
        return res

    def _excluded(self, A: Edge, B: Edge) -> set:
        key = (A.gid, B.gid)
        s = self._excl.get(key)
        if s is None:
            s = self._neighbors(A) | self._neighbors(B)
            self._excl[key] = s
        return s

    def probe(self, A: Edge, t: float):
        """Луч из A(t) внутрь материала. Возвращает (B, tB, a, b, толщина) или None."""
        o = A.point(t)
        d = A.normal(t)
        best = None
        for e in self.edges:
            if e is A or not _ray_may_hit(o, d, e.bbox):
                continue
            h = _ray_edge(e, o, d, self.eps)
            if h is not None and (best is None or h[0] < best[0]):
                best = (h[0], e, h[1])
        if best is None:
            return None
        _, B, tb = best
        if B.kind != A.kind and not self.prm.allow_mixed:
            return None  # прямая с дугой: ошибочный штрих, оси не делаем
        nB = B.normal(tb)
        if -(d[0] * nB[0] + d[1] * nB[1]) < self.cos_dev:
            return None  # стороны не обращены друг к другу
        b = B.point(tb)
        w = dist(o, b)
        if w < 2 * self.eps:
            return None
        # A и B замыкаются общими соседями в четырёхугольник, и обе короче расстояния между ними:
        # это два торца короткой стойки (i, l, перекладина), а не стороны штриха
        if max(A.length, B.length) < 0.9 * w and len((self._neighbors(A) & self._neighbors(B)) - {A.gid, B.gid}) >= 2:
            return None
        m = ((o[0] + b[0]) / 2, (o[1] + b[1]) / 2)
        thr = self.prm.disc_slack * w / 2.0
        thr2 = thr * thr
        excl = self._excluded(A, B)
        for e in self.edges:  # круг с диаметром o-b должен быть пуст
            if e.gid in excl:
                continue
            bb = e.bbox
            dx = max(bb[0] - m[0], 0.0, m[0] - bb[2])
            dy = max(bb[1] - m[1], 0.0, m[1] - bb[3])
            if dx * dx + dy * dy >= thr2:
                continue
            if dist_to_polyline(m, e.poly) < thr:
                return None
        return (B, tb, o, b, w)

    @staticmethod
    def _key(s) -> int:
        return -1 if s is None else s[0].gid

    def scan_edge(self, A: Edge) -> List[list]:
        """Серии подряд идущих лучей из A, попавших в одно и то же ребро."""
        step = self.w_est / (5.0 if A.kind == "line" else 6.0)  # шаг от толщины, а не от размера em
        n = max(6, min(60, int(math.ceil(A.length / step))))
        ts = [i / n for i in range(n + 1)]
        st = [self.probe(A, t) for t in ts]
        seq = []
        for i in range(n + 1):
            seq.append((ts[i], st[i]))
            if i < n and self._key(st[i]) != self._key(st[i + 1]):
                lo, slo, hi, shi = ts[i], st[i], ts[i + 1], st[i + 1]
                klo = self._key(slo)
                for _ in range(8):  # уточняем границу делением пополам
                    tm = (lo + hi) / 2
                    sm = self.probe(A, tm)
                    if self._key(sm) == klo:
                        lo, slo = tm, sm
                    else:
                        hi, shi = tm, sm
                seq.append((lo, slo))
                seq.append((hi, shi))
        runs, cur = [], []
        for t, s in seq:
            k = self._key(s)
            if k < 0 or (cur and self._key(cur[-1][1]) != k):
                if len(cur) > 1:
                    runs.append(cur)
                cur = []
            if k >= 0:
                cur.append((t, s))
        if len(cur) > 1:
            runs.append(cur)
        return runs

    # --- куски -------------------------------------------------------------------
    def make_piece(self, A: Edge, run: list) -> Optional[Piece]:
        B = run[0][1][0]
        samples: List[Sample] = [(t, s[1], s[2], s[3], s[4]) for t, s in run]
        m = [((s[2][0] + s[3][0]) / 2, (s[2][1] + s[3][1]) / 2) for s in samples]
        W = sum(s[4] for s in samples) / len(samples)
        La = sum(dist(m[i], m[i + 1]) for i in range(len(m) - 1))
        kind, axis = fit_axis(m, self.prm.flat_tol)
        axis_parts = ([axis] if kind == "line" else
                      fit_axis_segments(m, self.prm.flat_tol, max(1.5 * self.prm.flat_tol, 0.05 * W)))
        # короткий кусок обычно пятно (середины лучей круглой точки сливаются, направление оси
        # неизвестно). Но у пары почти параллельных прямых сторон направление оси известно точно,
        # поэтому зажатая между стойками диагональ (и, й, N) остаётся штрихом
        parallel = A.kind == "line" and B.kind == "line" and abs(vcross(A.t_in, B.t_in)) < 0.26
        if La < (0.1 if parallel else 0.3) * W:
            kind = "blob"
        extent = max(A.length * (samples[-1][0] - samples[0][0]),
                     B.length * (max(s[1] for s in samples) - min(s[1] for s in samples)))
        if kind == "blob" and extent < 0.2 * W:
            return None  # шум: лучи у самой границы рёбер
        t0 = end_tangents(axis_parts[0])[0]
        t1 = end_tangents(axis_parts[-1])[1]
        tbs = [s[1] for s in samples]
        return Piece(A, B, samples, kind, axis, W, La, t0, t1,
                     (samples[0][0], samples[-1][0]), (min(tbs), max(tbs)), axis_parts)

    def collect(self) -> List[Piece]:
        pieces: List[Piece] = []
        for A in self.edges:
            for run in self.scan_edge(A):
                p = self.make_piece(A, run)
                if p is not None:
                    pieces.append(p)
        # (A->B) и (B->A) описывают один и тот же кусок: оставляем длинный
        order = sorted(range(len(pieces)), key=lambda i: -pieces[i].length)
        keep = [True] * len(pieces)
        for ii, i in enumerate(order):
            if not keep[i]:
                continue
            P = pieces[i]
            for j in order[ii + 1:]:
                Q = pieces[j]
                if keep[j] and Q.a is P.b and Q.b is P.a:
                    lo, hi = max(P.ta[0], Q.tb[0]), min(P.ta[1], Q.tb[1])
                    short = min(P.ta[1] - P.ta[0], Q.tb[1] - Q.tb[0])
                    if hi - lo >= 0.5 * short:
                        keep[j] = False
        return [p for p, k in zip(pieces, keep) if k]

    # --- склейка -----------------------------------------------------------------
    @staticmethod
    def _outward(P: Piece, which: int) -> Pt:
        return (-P.t_start[0], -P.t_start[1]) if which == 0 else P.t_end

    @staticmethod
    def _turn(P: Piece) -> int:
        c = vcross(P.t_start, P.t_end)
        return 0 if abs(c) < 0.02 else (1 if c > 0 else -1)

    def _gap_inside(self, p: Pt, q: Pt) -> bool:
        """Проверяет весь предполагаемый мост, а не только три фиксированные точки.

        Шаг выборки привязан к оценённой толщине штриха: узкая полость или вырез
        между тремя прежними контрольными точками больше не остаётся незамеченным.
        """
        length = dist(p, q)
        step = max(1.0, 0.25 * self.w_est)
        n = max(3, min(64, int(math.ceil(length / step))))
        for i in range(1, n):
            if winding_number(lerp(p, q, i / n), self.a.polys) == 0:
                return False
        return True

    def compatible(self, P: Piece, Q: Piece, allow_stub: bool = False) -> bool:
        """allow_stub=False: склеиваются только прямая с прямой и дуга с дугой.
        allow_stub=True: допускается дуга с короткой прямой вставкой."""
        prm = self.prm
        wmax = max(P.width, Q.width)
        best = None
        for kp in (0, 1):
            pe = P.axis[0] if kp == 0 else P.axis[3]
            for kq in (0, 1):
                qe = Q.axis[0] if kq == 0 else Q.axis[3]
                d = dist(pe, qe)
                if best is None or d < best[0]:
                    best = (d, kp, kq, pe, qe)
        d, kp, kq, pe, qe = best
        # толщина должна быть непрерывной именно в месте стыка (у контрастных дуг
        # она вдоль штриха меняется в разы, поэтому средние значения не сравниваем)
        wp = P.samples[0][4] if kp == 0 else P.samples[-1][4]
        wq = Q.samples[0][4] if kq == 0 else Q.samples[-1][4]
        if max(wp, wq) > 1.5 * min(wp, wq):
            return False
        if P.kind == "line" and Q.kind == "line" and not allow_stub:
            u = vunit(vsub(P.axis[3], P.axis[0]))
            v = vunit(vsub(Q.axis[3], Q.axis[0]))
            if not u or not v or abs(vcross(u, v)) > math.sin(math.radians(prm.line_merge_deg)):
                return False
            uu = u if P.length >= Q.length else v
            g = vsub(qe, pe)
            if abs(vdot(g, uu)) > prm.gap_factor * wmax:
                return False
            ends4 = [P.axis[0], P.axis[3], Q.axis[0], Q.axis[3]]
            f0, f1 = max(((x, y) for x in ends4 for y in ends4), key=lambda xy: dist(xy[0], xy[1]))
            w = vunit(vsub(f1, f0))
            if not w or any(abs(vcross(vsub(q, f0), w)) > 0.5 * wmax for q in ends4):
                return False  # объединение не ложится на одну прямую (стороны V, A, M)
        else:  # дуга с дугой или (allow_stub) дуга с короткой прямой вставкой
            if (P.kind == "line") != (Q.kind == "line") and not allow_stub:
                return False
            if P.kind == "line" and Q.kind == "line":
                return False
            if d > wmax:
                return False
            if abs(angle_deg(self._outward(P, kp), (-self._outward(Q, kq)[0], -self._outward(Q, kq)[1]))) > prm.merge_deg:
                return False
            if P.kind == "arc" and Q.kind == "arc":
                sp = self._turn(P) * (1 if kp == 1 else -1)
                sq = self._turn(Q) * (1 if kq == 0 else -1)
                if sp * sq < 0:
                    return False
        if d > 0.25 * wmax and not self._gap_inside(pe, qe):
            return False
        return True

    @staticmethod
    def _chains(pieces: List[Piece], w: float, merge_deg: float = 45.0) -> List[Tuple[List[Piece], List[Cubic], List[Piece]]]:
        """Собирает только действительно примыкающие фрагменты одной дуги.

        Расстояние само по себе не считается достаточным основанием для склейки:
        дополнительно проверяются направление касательных и локальная толщина.
        Разрыв, который нельзя безопасно закрыть, остаётся отдельной осью и может
        быть обработан позднее общим соединителем.
        """
        rest = list(pieces)
        out: List[Tuple[List[Piece], List[Cubic], List[Piece]]] = []
        max_gap = 0.35 * max(w, 1e-6)
        angle_limit = max(1.0, merge_deg)

        def ends(p):
            return (p.axis[0], p.axis[3])

        def oriented_parts(p: Piece, reverse: bool = False) -> List[Cubic]:
            parts = p.axis_parts or [p.axis]
            return [_reverse_cubic(c) for c in reversed(parts)] if reverse else list(parts)

        def endpoint_width(p: Piece, at_start: bool) -> float:
            return p.samples[0][4] if at_start else p.samples[-1][4]

        while rest:
            if len(rest) == 1:
                out.append(([rest[0]], [rest[0].axis], []))
                break
            best = None
            for p in rest:  # начало цепочки: самый свободный конец
                for k, e in enumerate(ends(p)):
                    dm = min(dist(e, q) for r in rest if r is not p for q in ends(r))
                    if best is None or dm > best[0]:
                        best = (dm, p, k)
            _, start, k = best
            first_parts = oriented_parts(start, reverse=(k == 1))
            chain_p = [start]
            chain_c = list(first_parts)
            chain_end_width = endpoint_width(start, at_start=(k == 1))
            rest.remove(start)

            while True:
                E = chain_c[-1][3]
                tangent_end = end_tangents(chain_c[-1])[1]
                options = []
                for r in rest:
                    for kk, e in enumerate(ends(r)):
                        gap = dist(E, e)
                        if gap > max_gap:
                            continue
                        next_parts = oriented_parts(r, reverse=(kk == 1))
                        tangent_next = end_tangents(next_parts[0])[0]
                        turn = abs(angle_deg(tangent_end, tangent_next))
                        if turn > angle_limit:
                            continue
                        join_width = endpoint_width(r, at_start=(kk == 0))
                        far_width = endpoint_width(r, at_start=(kk == 1))
                        if min(chain_end_width, join_width) <= 1e-6:
                            continue
                        width_ratio = max(chain_end_width, join_width) / min(chain_end_width, join_width)
                        if width_ratio > 1.5:
                            continue
                        # Выбираем не просто ближайший фрагмент, а наиболее
                        # согласованный по расстоянию, касательной и толщине.
                        score = (gap / max(w, 1e-6)
                                 + 0.25 * turn / angle_limit
                                 + 0.15 * abs(math.log(width_ratio)))
                        options.append((score, r, kk, next_parts, far_width))
                if not options:
                    break
                _, r, kk, next_parts, far_width = min(options, key=lambda item: item[0])
                chain_c.extend(next_parts)
                chain_p.append(r)
                chain_end_width = far_width
                rest.remove(r)

            poly = [q for c in chain_c for q in flatten(c, 8)]
            dups = [r for r in rest if all(dist_to_polyline(e, poly) <= 0.7 * w for e in ends(r))]
            for r in dups:
                rest.remove(r)
            out.append((chain_p, chain_c, dups))
        return out

    @staticmethod
    def _patch(p: Piece) -> List[Pt]:
        return [s[2] for s in p.samples] + [s[3] for s in reversed(p.samples)]

    def build(self, pieces: List[Piece]) -> List[Stroke]:
        prm = self.prm
        cand = [p for p in pieces if p.kind != "blob"]
        parent = list(range(len(cand)))

        def find(i):
            while parent[i] != i:
                parent[i] = parent[parent[i]]
                i = parent[i]
            return i

        def collinear(idxs: List[int]) -> bool:
            """Концы осей всех кусков группы лежат около одной прямой."""
            pts = [e for k in idxs for e in (cand[k].axis[0], cand[k].axis[3])]
            f0, f1 = max(((x, y) for x in pts for y in pts), key=lambda xy: dist(xy[0], xy[1]))
            w = vunit(vsub(f1, f0))
            if not w:
                return True
            wm = sum(cand[k].width for k in idxs) / len(idxs)
            if not all(abs(vcross(vsub(q, f0), w)) <= 0.6 * wm for q in pts):
                return False
            # и ось каждого достаточно длинного куска смотрит в ту же сторону, что общая прямая
            # (на большом зазоре отклонение в 0.6 толщины почти не видно, а кусок уже повёрнут)
            lim = math.sin(math.radians(1.5 * self.prm.line_merge_deg))
            for k in idxs:
                dk = vunit(vsub(cand[k].axis[3], cand[k].axis[0]))
                if dk and cand[k].length >= 0.5 * cand[k].width and abs(vcross(dk, w)) > lim:
                    return False
            return True

        # проход 1: прямая с прямой, дуга с дугой. Пары по возрастанию зазора; объединение прямых
        # принимается, только если вся группа остаётся на одной прямой (иначе цепочка кусков,
        # каждый из которых чуть повёрнут к соседу, уходит по диагонали: перекладина f с завитком)
        pairs = []
        for i in range(len(cand)):
            for j in range(i + 1, len(cand)):
                if self.compatible(cand[i], cand[j]):
                    gap = min(dist(x, y) for x in (cand[i].axis[0], cand[i].axis[3])
                              for y in (cand[j].axis[0], cand[j].axis[3]))
                    pairs.append((gap, i, j))
        pairs.sort()
        members: Dict[int, List[int]] = {i: [i] for i in range(len(cand))}
        for _, i, j in pairs:
            ri, rj = find(i), find(j)
            if ri == rj:
                continue
            merged = members[ri] + members[rj]
            if all(cand[k].kind == "line" for k in merged) and not collinear(merged):
                continue
            parent[rj] = ri
            members[ri] = merged
            del members[rj]

        def group_stats(idxs):
            ps = [cand[i] for i in idxs]
            ws = [s[4] for p in ps for s in p.samples]
            if all(p.kind != "arc" for p in ps):
                pts = [e for p in ps for e in (p.axis[0], p.axis[3])]
                ln = max(dist(x, y) for x in pts for y in pts)
            else:
                ln = sum(p.length for p in ps)
            return ln, sum(ws) / len(ws)

        # проход 2: совсем короткая прямая между дугами (спинка S) входит в дугу
        by_root: Dict[int, List[int]] = {}
        for i in range(len(cand)):
            by_root.setdefault(find(i), []).append(i)
        for root, idxs in list(by_root.items()):
            if any(cand[i].kind != "line" for i in idxs):
                continue
            ln, wm = group_stats(idxs)
            if ln > 1.0 * wm:
                continue
            best = None
            for i in idxs:
                for j in range(len(cand)):
                    if cand[j].kind == "line" or find(j) == find(i):
                        continue
                    if self.compatible(cand[i], cand[j], allow_stub=True):
                        gap = min(dist(x, y) for x in (cand[i].axis[0], cand[i].axis[3])
                                  for y in (cand[j].axis[0], cand[j].axis[3]))
                        if best is None or gap < best[0]:
                            best = (gap, j)
            if best is not None:
                parent[find(idxs[0])] = find(best[1])
        groups: Dict[int, List[Piece]] = {}
        for i, p in enumerate(cand):
            groups.setdefault(find(i), []).append(p)

        # микрокуски (очень короткие рёбра у экстремумов контрастных дуг) присоединяем
        # к ближайшему концу штриха, если толщина в месте стыка непрерывна
        loose: List[Piece] = []
        for p in pieces:
            if p.kind != "blob":
                continue
            best = None
            for i, q in enumerate(cand):
                for kp in (0, 1):
                    pe = p.axis[0] if kp == 0 else p.axis[3]
                    wp = p.samples[0][4] if kp == 0 else p.samples[-1][4]
                    for kq in (0, 1):
                        qe = q.axis[0] if kq == 0 else q.axis[3]
                        wq = q.samples[0][4] if kq == 0 else q.samples[-1][4]
                        d = dist(pe, qe)
                        if d <= 0.3 * max(wp, wq) and max(wp, wq) <= 1.5 * min(wp, wq) \
                                and (best is None or d < best[0]):
                            best = (d, i)
            if best is None:
                loose.append(p)
            else:
                groups[find(best[1])].append(p)
        pieces = [p for p in pieces if p.kind != "blob"] + loose

        strokes: List[Stroke] = []
        blobs: List[Tuple[List[Piece], List[Pt]]] = []  # (куски, точки)
        for ps in groups.values():
            ws = [s[4] for p in ps for s in p.samples]
            w_min, w_max, w_mean = min(ws), max(ws), sum(ws) / len(ws)
            if all(p.kind != "arc" for p in ps):
                kind = "line"
                pts = [e for p in ps if p.kind != "blob" for e in (p.axis[0], p.axis[3])]
                a0, a1 = max(((x, y) for x in pts for y in pts), key=lambda xy: dist(xy[0], xy[1]))
                axis = [line_cubic(a0, a1)]
                length = dist(a0, a1)
            else:
                blob_ps = [p for p in ps if p.kind == "blob"]
                for ci, (cp, cc, dups) in enumerate(
                        self._chains([p for p in ps if p.kind != "blob"], w_mean, prm.merge_deg)):
                    pcs = cp + dups + (blob_ps if ci == 0 else [])
                    cws = [q[4] for p in pcs for q in p.samples]
                    cmean = sum(cws) / len(cws)
                    clen = sum(_cubic_len(c) for c in cc)
                    if prm.blob_ratio > 0 and clen < prm.blob_ratio * cmean:
                        blobs.append((pcs, [q for p in pcs for sm in p.samples for q in (sm[2], sm[3])]))
                        continue
                    join_tolerance = 0.35 * cmean
                    closed = len(cc) > 1 and dist(cc[0][0], cc[-1][3]) <= join_tolerance
                    axis = close_gaps(cc, closed, join_tolerance)
                    if prm.smooth_joints:
                        axis = smooth_chain(axis, closed, prm.snap_deg)
                    strokes.append(Stroke(0, "arc", pcs, axis, [self._patch(p) for p in pcs],
                                          min(cws), max(cws), cmean, clen, closed))
                continue
            if prm.blob_ratio > 0 and length < prm.blob_ratio * w_mean:
                blobs.append((ps, [q for p in ps for s in p.samples for q in (s[2], s[3])]))
                continue
            strokes.append(Stroke(0, kind, ps, axis, [self._patch(p) for p in ps],
                                  w_min, w_max, w_mean, length, False))
        for p in pieces:
            if p.kind == "blob":
                blobs.append(([p], [q for s in p.samples for q in (s[2], s[3])]))

        # близкие пятна объединяем
        def bb(pts):
            xs, ys = [q[0] for q in pts], [q[1] for q in pts]
            return min(xs) - 1, min(ys) - 1, max(xs) + 1, max(ys) + 1

        bparent = list(range(len(blobs)))

        def bfind(i):
            while bparent[i] != i:
                bparent[i] = bparent[bparent[i]]
                i = bparent[i]
            return i

        boxes = [bb(b[1]) for b in blobs]
        for i in range(len(blobs)):
            for j in range(i + 1, len(blobs)):
                x0, y0 = max(boxes[i][0], boxes[j][0]), max(boxes[i][1], boxes[j][1])
                x1, y1 = min(boxes[i][2], boxes[j][2]), min(boxes[i][3], boxes[j][3])
                if x1 > x0 and y1 > y0:
                    ov = (x1 - x0) * (y1 - y0)
                    ar = lambda b: (b[2] - b[0]) * (b[3] - b[1])
                    if ov >= 0.3 * min(ar(boxes[i]), ar(boxes[j])):
                        bparent[bfind(j)] = bfind(i)
        bgroups: Dict[int, List[int]] = {}
        for i in range(len(blobs)):
            bgroups.setdefault(bfind(i), []).append(i)
        for idxs in bgroups.values():
            ps = [p for i in idxs for p in blobs[i][0]]
            pts = [q for i in idxs for q in blobs[i][1]]
            ws = [s[4] for p in ps for s in p.samples]
            strokes.append(Stroke(0, "blob", ps, [], [convex_hull(pts)],
                                  min(ws), max(ws), sum(ws) / len(ws), 0.0))

        # пятно, лежащее почти целиком внутри настоящего штриха (пара торцов
        # прямоугольника), только мешает: убираем
        solid = [patch for s in strokes if s.kind != "blob" for patch in s.patches if len(patch) >= 3]
        if solid:
            kept = []
            for s in strokes:
                if s.kind == "blob" and len(s.patches[0]) >= 3:
                    hull = s.patches[0]
                    cxy = (sum(q[0] for q in hull) / len(hull), sum(q[1] for q in hull) / len(hull))
                    probes = [lerp(cxy, q, 0.85) for q in hull] + [cxy]
                    inside = sum(1 for q in probes if any(winding_number(q, [pa]) != 0 for pa in solid))
                    if inside >= 0.75 * len(probes):
                        continue
                kept.append(s)
            strokes = kept

        def cx(s: Stroke) -> float:
            pts = [q for patch in s.patches for q in patch]
            return sum(q[0] for q in pts) / len(pts)

        strokes.sort(key=cx)
        for i, s in enumerate(strokes, 1):
            s.sid = i
            s.hue = ((i * 137.508) % 360.0) / 360.0
            s.tags = ["fit"] * len(s.axis)
            s.prof = [(((q[2][0] + q[3][0]) / 2, (q[2][1] + q[3][1]) / 2), q[4])
                      for p in s.pieces for q in p.samples]
        return strokes


def mark_dots(a: "Analysis", prm: Params) -> None:
    """Отдельная компактная деталь (точка над i и j, двоеточие, период), размер которой не больше
    dot_ratio толщин, это не набор штрихов, а точка. Штрихи внутри неё заменяются одним
    объектом вида 'dot'; тонкое начертание получает ту же точку, уменьшенную вокруг центра."""
    if prm.dot_ratio <= 0:
        return
    n = len(a.strokes)
    parent = list(range(n))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    owner: Dict[int, int] = {}
    for i, s in enumerate(a.strokes):
        for p in s.pieces:
            for c in (p.a.contour, p.b.contour):
                if c in owner:
                    parent[find(i)] = find(owner[c])
                else:
                    owner[c] = i
    groups: Dict[int, List[int]] = {}
    for i in range(n):
        groups.setdefault(find(i), []).append(i)
    real = [s.w_mean for s in a.strokes if s.kind in ("line", "arc")]
    w_ref = sorted(real)[len(real) // 2] if real else None  # типичная толщина штриха этой буквы
    out: List[Stroke] = []
    for idxs in sorted(groups.values(), key=lambda g: g[0]):
        ss = [a.strokes[i] for i in idxs]
        conts = {c for s in ss for p in s.pieces for c in (p.a.contour, p.b.contour)}
        wm = sum(s.w_mean for s in ss) / len(ss)
        if len(conts) == 1:
            c = next(iter(conts))
            xs = [q[0] for q in a.polys[c]]
            ys = [q[1] for q in a.polys[c]]
            if max(max(xs) - min(xs), max(ys) - min(ys)) <= prm.dot_ratio * wm:
                d = Stroke(0, "dot", [p for s in ss for p in s.pieces], [], [list(a.polys[c])],
                           min(s.w_min for s in ss), max(s.w_max for s in ss), wm, 0.0)
                d.contour = c
                out.append(d)
                continue
        out.extend(ss)
    # контур, не принадлежащий ни одному штриху (круглая точка не даёт пар сторон), но компактный
    for c, poly in enumerate(a.polys):
        if c in owner or len(poly) < 3:
            continue
        xs = [q[0] for q in poly]
        ys = [q[1] for q in poly]
        w, h = max(xs) - min(xs), max(ys) - min(ys)
        if min(w, h) <= 0 or max(w, h) / min(w, h) > 2.0:
            continue
        if w_ref is not None and max(w, h) > prm.dot_ratio * w_ref:
            continue
        d = Stroke(0, "dot", [], [], [list(poly)], min(w, h), max(w, h), (w + h) / 2, 0.0)
        d.contour = c
        out.append(d)

    def cx(s: Stroke) -> float:
        pts = [q for patch in s.patches for q in patch]
        return sum(q[0] for q in pts) / len(pts)

    out.sort(key=cx)
    for i, s in enumerate(out, 1):
        s.sid = i
        s.hue = ((i * 137.508) % 360.0) / 360.0
    a.strokes = out


def find_strokes(a: "Analysis", upm: int, prm: Params) -> None:
    """Заполняет a.strokes и покрытие рёбер (Edge.cover)."""
    for e in a.edges:
        e.cover = []
    sf = _StrokeFinder(a, upm, prm)
    pieces = sf.collect()
    for p in pieces:
        p.a.cover.append(p.ta)
        p.b.cover.append(p.tb)
    a.strokes = sf.build(pieces)
    mark_dots(a, prm)


# ============================================================================
#  Часть 3в. Общая ось: соединение осей штрихов, предсказание контура
# ============================================================================
#
#  Опора на существующий контур (чтобы «то, что есть, имело вес»):
#   1. Конец штриха: луч вперёд сразу упирается в контур, значит, это настоящий торец,
#      и соединять его не с чем. Если вперёд идёт материал, это стык.
#   2. Соединение принимается только при опоре в контуре: рядом с точкой стыка должны
#      быть вогнутые углы контура (как у T, L, X, у стойки с аркой).
#   3. Продление оси идёт только внутри материала. Свободный конец, перед которым
#      виден торец, достраивается до этого торца.
#   4. Проверка: по осям и измеренной толщине строится «лента» штриха. При исходной
#      толщине она должна повторять существующий контур, а при тонкой даёт
#      предсказание тонкого начертания. Расхождение с контуром выводится числами.


@dataclass(eq=False)
class Joint:
    kind: str  # 'T' | 'L' | 'Y' | 'X' | 'G1'
    pos: Pt
    strokes: List["Stroke"]
    support: int  # число вогнутых углов контура рядом со стыком


@dataclass(eq=False)
class End:
    s: "Stroke"
    idx: int  # 0: начало оси, 1: конец оси
    pos: Pt
    t: Pt  # единичная касательная наружу
    w: float  # толщина штриха у конца
    terminal: bool = False
    used: bool = False
    mark: str = ""  # 'terminal' | 'cap' | 'open' | ''


def _cubic_poly(c: Cubic) -> List[Pt]:
    return [c[0], c[3]] if is_flat(c, 0.05) else flatten(c, 12)


def _nearest_on_segment(p: Pt, a: Pt, b: Pt) -> Pt:
    ex, ey = b[0] - a[0], b[1] - a[1]
    l2 = ex * ex + ey * ey
    t = 0.0 if l2 == 0 else max(0.0, min(1.0, ((p[0] - a[0]) * ex + (p[1] - a[1]) * ey) / l2))
    return (a[0] + ex * t, a[1] + ey * t)


def _axis_poly(s: "Stroke") -> List[Pt]:
    pts: List[Pt] = []
    for c in s.axis:
        f = _cubic_poly(c)
        if pts and dist(pts[-1], f[0]) < 1e-6:
            f = f[1:]
        pts.extend(f)
    return pts


def _axis_samples(s: "Stroke", step: float) -> List[Pt]:
    pts: List[Pt] = []
    for c in s.axis:
        n = max(2, int(math.ceil(_cubic_len(c) / step)))
        if not is_flat(c, 0.05):
            n = max(n, 8)
        f = flatten(c, n)
        if pts and dist(pts[-1], f[0]) < 1e-6:
            f = f[1:]
        pts.extend(f)
    return pts


def _width_at(s: "Stroke", pt: Pt) -> float:
    """Толщина штриха в точке оси: по двум ближайшим измеренным отрезкам A-B."""
    prof = s.prof
    if not prof:
        return s.w_mean
    best = sorted(prof, key=lambda m: (m[0][0] - pt[0]) ** 2 + (m[0][1] - pt[1]) ** 2)[:2]
    if len(best) == 1:
        return best[0][1]
    d0, d1 = dist(best[0][0], pt), dist(best[1][0], pt)
    if d0 + d1 < 1e-9:
        return best[0][1]
    return (best[0][1] * d1 + best[1][1] * d0) / (d0 + d1)


def _seg_ray(p: Pt, d: Pt, a: Pt, b: Pt) -> Optional[Tuple[float, float]]:
    """Луч p + s*d и отрезок a-b: (s, u) или None."""
    ex, ey = b[0] - a[0], b[1] - a[1]
    den = d[0] * ey - d[1] * ex
    if abs(den) < 1e-12:
        return None
    wx, wy = a[0] - p[0], a[1] - p[1]
    return (wx * ey - wy * ex) / den, (wx * d[1] - wy * d[0]) / den


def _orient(poly: List[Pt], ccw: bool) -> List[Pt]:
    area = sum(poly[i][0] * poly[(i + 1) % len(poly)][1] - poly[(i + 1) % len(poly)][0] * poly[i][1]
               for i in range(len(poly))) / 2.0
    return poly if (area >= 0) == ccw else poly[::-1]


def _resample(poly: List[Pt], step: float) -> List[Pt]:
    """Точки через равные промежутки вдоль замкнутой ломаной."""
    out: List[Pt] = []
    n = len(poly)
    acc = 0.0
    for i in range(n):
        a, b = poly[i], poly[(i + 1) % n]
        L = dist(a, b)
        pos = acc
        while pos < L:
            out.append(lerp(a, b, pos / L))
            pos += step
        acc = pos - L
    return out



# ----------------------------------------------------------------------------
#  Срединная ось по настоящему контуру.
#  Точка оси в стыке — центр круга, который касается контура с двух сторон. Расстояние до
#  контура считается по векторным рёбрам (прямые и кубики, с уточнением Ньютоном), так что
#  вогнутые углы стыка (вырез между стойкой и чашей) работают как отдельные точки-опоры.
#  Ось ведётся шагами: предсказание по касательной, затем поправка вдоль линии между
#  ближайшими точками контура слева и справа, пока расстояния не сравняются.
#  Растр и диаграмма Вороного не нужны.
# ----------------------------------------------------------------------------


def _bez_d2(p: Cubic, t: float) -> Pt:
    mt = 1.0 - t
    a = (p[2][0] - 2 * p[1][0] + p[0][0], p[2][1] - 2 * p[1][1] + p[0][1])
    b = (p[3][0] - 2 * p[2][0] + p[1][0], p[3][1] - 2 * p[2][1] + p[1][1])
    return (6.0 * (mt * a[0] + t * b[0]), 6.0 * (mt * a[1] + t * b[1]))


def _edge_nearest(e: "Edge", c: Pt) -> Tuple[float, Pt]:
    """Ближайшая к c точка ребра: (расстояние, точка). Для кривой уточняется по настоящей кубике."""
    poly = e.poly
    n = len(poly) - 1
    best_d, best_q, best_t = 1e18, poly[0], 0.0
    for i in range(n):
        a, b = poly[i], poly[i + 1]
        abx, aby = b[0] - a[0], b[1] - a[1]
        l2 = abx * abx + aby * aby
        u = 0.0 if l2 == 0 else max(0.0, min(1.0, ((c[0] - a[0]) * abx + (c[1] - a[1]) * aby) / l2))
        q = (a[0] + abx * u, a[1] + aby * u)
        d = math.hypot(c[0] - q[0], c[1] - q[1])
        if d < best_d:
            best_d, best_q, best_t = d, q, (i + u) / n
    if e.kind == "curve":
        t = best_t
        for _ in range(4):
            B = bez_pt(e.pts, t)
            D1 = bez_deriv(e.pts, t)
            D2 = _bez_d2(e.pts, t)
            f = (B[0] - c[0]) * D1[0] + (B[1] - c[1]) * D1[1]
            fp = D1[0] * D1[0] + D1[1] * D1[1] + (B[0] - c[0]) * D2[0] + (B[1] - c[1]) * D2[1]
            if abs(fp) < 1e-9:
                break
            t2 = min(1.0, max(0.0, t - f / fp))
            B2 = bez_pt(e.pts, t2)
            d2 = math.hypot(c[0] - B2[0], c[1] - B2[1])
            if d2 > best_d + 1e-9:
                break
            t, best_q, best_d = t2, B2, d2
    return best_d, best_q


def _fit_path(pts: List[Pt], t_start: Optional[Pt], t_end: Optional[Pt], tol: float) -> List[Cubic]:
    """Кубики Безье вдоль трассы. t_start, t_end: касательные вперёд по пути (None: по точкам)."""
    p: List[Pt] = [pts[0]]
    for q in pts[1:]:
        if dist(q, p[-1]) > 1e-6:
            p.append(q)
    if len(p) < 2:
        return []
    chord = vsub(p[-1], p[0])
    L = math.hypot(chord[0], chord[1])
    straight = len(p) == 2 or (
        L > 1e-9 and all(abs(vcross(chord, vsub(q, p[0]))) / L <= tol for q in p[1:-1]))
    if straight and L > 1e-9:
        u = (chord[0] / L, chord[1] / L)
        if (t_start is None or abs(vcross(t_start, u)) < 0.05) and (t_end is None or abs(vcross(t_end, u)) < 0.05):
            return [line_cubic(p[0], p[-1])]
    t1 = t_start or vunit(vsub(p[1], p[0])) or (1.0, 0.0)
    te = t_end or vunit(vsub(p[-1], p[-2])) or (1.0, 0.0)
    out: List[Cubic] = []
    _fit_cubic(p, t1, (-te[0], -te[1]), tol, out)
    return out


def _arc_cubics(C: Pt, rho: float, P: Pt, t: Pt, sigma: int, theta: float) -> List[Cubic]:
    """Дуга окружности (центр C, радиус rho) от точки P с касательной t, поворот на theta (sigma=+-1: сторона)."""
    n = max(1, int(math.ceil(theta / (math.pi / 2) - 1e-9)))
    d = theta / n
    k = 4.0 / 3.0 * math.tan(d / 4.0) * rho
    ca, sa = math.cos(sigma * d), math.sin(sigma * d)
    out: List[Cubic] = []
    p0, t0 = P, t
    for _ in range(n):
        r = (p0[0] - C[0], p0[1] - C[1])
        p3 = (C[0] + r[0] * ca - r[1] * sa, C[1] + r[0] * sa + r[1] * ca)
        t3 = (t0[0] * ca - t0[1] * sa, t0[0] * sa + t0[1] * ca)
        out.append((p0, (p0[0] + t0[0] * k, p0[1] + t0[1] * k), (p3[0] - t3[0] * k, p3[1] - t3[1] * k), p3))
        p0, t0 = p3, t3
    return out


@dataclass(eq=False)
class _Link:
    """Готовый участок оси, соединяющий конец штриха со стыком."""
    cubics: List[Cubic]
    tag: str = "ext"
    prof: List[Tuple[Pt, float]] = field(default_factory=list)  # (точка оси, толщина) по радиусу вписанного круга
    cut: Optional[Pt] = None  # перед присоединением ось штриха обрезается в этой точке (участок заменяется)
    how: str = ""  # способ: 'oval' (дуга чаши до стойки) или 'medial' (середина контура)

    @property
    def end(self) -> Pt:
        return self.cubics[-1][3]


class _Connector:
    G1_DEG = 20.0  # концы почти продолжают друг друга: плавное соединение

    def __init__(self, a: "Analysis", upm: int, prm: Params):
        self.a, self.upm, self.prm = a, upm, prm
        self.strokes = [s for s in a.strokes if s.axis]
        self.joints: List[Joint] = []
        self.corners: List[Tuple[Pt, End, End]] = []
        self.polys: Dict[int, List[Pt]] = {}
        self.ink = {id(s): {c for p in s.pieces for c in (p.a.contour, p.b.contour)} for s in self.strokes}
        self._memo: Dict[tuple, Optional[_Link]] = {}  # результаты трассировки: plan() зовётся несколько раз

    def record_cap(self, e: "End", h) -> None:
        """Торец, в который упирается конец оси: лента режется по его линии, а не перпендикулярно оси."""
        d = h[1].tangent(h[2])
        e.s.caps[e.idx] = (h[3], d)

    # --- опора в контуре ---------------------------------------------------------
    def inside(self, p: Pt) -> bool:
        return winding_number(p, self.a.polys) != 0

    def seg_inside(self, p: Pt, q: Pt, n: int = 4) -> bool:
        return all(self.inside(lerp(p, q, (i + 1) / (n + 1))) for i in range(n))

    def near_ink(self, X: Pt, tol: float) -> bool:
        """Точка стыка лежит в материале или не дальше tol от контура (вершина V за обрезанным основанием не годится)."""
        if self.inside(X):
            return True
        for e in self.a.edges:
            bb = e.bbox
            if bb[0] - tol <= X[0] <= bb[2] + tol and bb[1] - tol <= X[1] <= bb[3] + tol \
                    and dist_to_polyline(X, e.poly) <= tol:
                return True
        return False

    def support(self, X: Pt, w1: float, w2: float, sin_t: float = 1.0) -> int:
        """Опора в контуре: сколько вогнутых углов рядом со стыком. У острого стыка углы
        лежат дальше, поэтому радиус растёт как 1/sin угла между осями."""
        R = 1.0 * math.hypot(w1, w2) / max(sin_t, 0.4)
        return sum(1 for v in self.a.vertices if v.kind == "concave" and dist(v.pos, X) <= R)

    def same_ink(self, s1: "Stroke", s2: "Stroke") -> bool:
        """Штрихи лежат на одном и том же контуре, то есть в одной связной части буквы."""
        return bool(self.ink[id(s1)] & self.ink[id(s2)])

    def ahead(self, p: Pt, t: Pt):
        """Первое пересечение луча из p вперёд с контуром: (расстояние, ребро, t, точка)."""
        o = (p[0] - t[0], p[1] - t[1])  # отступаем на 1 ед. назад, чтобы быть внутри материала
        best = None
        for e in self.a.edges:
            if not _ray_may_hit(o, t, e.bbox):
                continue
            h = _ray_edge(e, o, t, 0.05)
            if h is not None and (best is None or h[0] < best[0]):
                best = (h[0], e, h[1])
        if best is None:
            return None
        s, e, tb = best
        return s - 1.0, e, tb, e.point(tb)

    def hits(self, o: Pt, d: Pt, limit: float) -> List[float]:
        """Расстояния до последовательных пересечений луча с контуром (до limit)."""
        res: List[float] = []
        cur, total = o, 0.0
        for _ in range(6):
            best = None
            for e in self.a.edges:
                if not _ray_may_hit(cur, d, e.bbox):
                    continue
                h = _ray_edge(e, cur, d, 0.05)
                if h is not None and (best is None or h[0] < best):
                    best = h[0]
            if best is None or total + best > limit:
                break
            total += best
            res.append(total)
            cur = (cur[0] + d[0] * (best + 0.1), cur[1] + d[1] * (best + 0.1))
            total += 0.1
        return res

    def corridor_center(self, q: Pt, n: Pt, reach: float) -> Optional[float]:
        """Середина ближайшего к q куска чернил на прямой q + u*n (u от -reach до reach)."""
        bounds = sorted([-h for h in self.hits(q, (-n[0], -n[1]), reach)] + self.hits(q, n, reach))
        edges = [-reach] + bounds + [reach]
        state = self.inside(q)
        # состояние внутри промежутка, содержащего u=0
        k = sum(1 for b in bounds if b < 0)
        best = None
        for i in range(len(edges) - 1):
            ink = state if (i == k) else (state == ((i - k) % 2 == 0))
            if not ink:
                continue
            lo, hi = edges[i], edges[i + 1]
            gap = 0.0 if lo <= 0.0 <= hi else min(abs(lo), abs(hi))
            if hi - lo > 1.0 and (best is None or gap < best[0]):
                best = (gap, (lo + hi) / 2.0)
        return None if best is None else best[1]

    def corridor_path(self, e: End, X: Pt, target: "Stroke") -> Optional[Cubic]:
        """Путь от конца e к оси target по середине прохода между контурами (шейка у выреза)."""
        d = vunit(vsub(X, e.pos))
        if d is None:
            return None
        n = (-d[1], d[0])
        pts: List[Pt] = [e.pos]
        D = dist(e.pos, X)
        f_edge = 1.0 - 0.5 * _width_at(target, X) / max(D, 1e-9)  # дальше начинается тело самой цели
        fs = [f for f in (0.15, 0.3, 0.45, 0.6, 0.75, 0.9) if f <= f_edge - 0.1] or [0.5]
        for f in fs:
            q = lerp(e.pos, X, f)
            u = self.corridor_center(q, n, 1.2 * e.w)
            if u is None:
                return None
            pts.append((q[0] + n[0] * u, q[1] + n[1] * u))
        poly = self.polys[id(target)]
        end = min((_nearest_on_segment(pts[-1], a_, b_) for a_, b_ in zip(poly, poly[1:])),
                  key=lambda c: dist(c, pts[-1]))
        pts.append(end)
        # гладкая кривая: касательная в начале по концу штриха, в конце по направлению подхода;
        # по точкам центра прохода она только проверяется (МНК давал волну у стойки)
        t_end = vunit(vsub(end, pts[-2])) or e.t
        h = 0.4 * dist(e.pos, end)
        c: Cubic = (e.pos, (e.pos[0] + e.t[0] * h, e.pos[1] + e.t[1] * h),
                    (end[0] - t_end[0] * h, end[1] - t_end[1] * h), end)
        if all(self.inside(q) for q in flatten(c, 10)[1:-1]):
            return c
        kind, c = fit_axis(pts, self.prm.flat_tol)  # запасной вариант: подгонка по точкам прохода
        h = max(dist(c[0], c[1]), dist(c[0], c[3]) / 3.0)
        c = (c[0], (c[0][0] + e.t[0] * h, c[0][1] + e.t[1] * h), c[2], c[3])
        if not all(self.inside(q) for q in flatten(c, 10)[1:-1]):
            return None
        return c


    # --- срединная ось по настоящему контуру -----------------------------------------
    def _sites(self, c: Pt, cap: float):
        """Ближайшие к c точки контура по одной на ребро: [(расстояние, точка, ребро)] по возрастанию."""
        res = []
        for e in self.a.edges:
            bb = e.bbox
            dx = max(bb[0] - c[0], 0.0, c[0] - bb[2])
            dy = max(bb[1] - c[1], 0.0, c[1] - bb[3])
            if dx * dx + dy * dy > cap * cap:
                continue
            d, q = _edge_nearest(e, c)
            if d <= cap:
                res.append((d, q, e))
        res.sort(key=lambda s: s[0])
        return res

    def _snap(self, c: Pt, t: Pt, cap: float, iters: int = 8):
        """Сдвигает c вдоль линии «правая точка контура - левая точка контура», пока расстояния
        до контура слева и справа от направления t не сравняются: (c, dl, dr, sites) или None."""
        for _ in range(iters):
            S = self._sites(c, cap)
            left = right = None
            for d, q, _e in S:
                cr = vcross(t, vsub(q, c))
                if abs(cr) < 0.25 * d:  # точка впереди или позади, а не сбоку
                    continue
                if cr > 0:
                    if left is None:
                        left = (d, q)
                elif right is None:
                    right = (d, q)
                if left is not None and right is not None:
                    break
            if left is None or right is None:
                return None
            dl, ql = left
            dr, qr = right
            if abs(dl - dr) < 0.03:
                return c, dl, dr, S
            u = vunit(vsub(ql, qr))
            if u is None:
                return None
            den = vdot(u, vsub(ql, c)) / dl + vdot(u, vsub(c, qr)) / dr
            if den < 0.15:
                return None
            step = max(-0.35 * cap, min(0.35 * cap, (dl - dr) / den))
            c = (c[0] + u[0] * step, c[1] + u[1] * step)
        S = self._sites(c, cap)
        left = right = None
        for d, q, _e in S:
            cr = vcross(t, vsub(q, c))
            if abs(cr) < 0.25 * d:
                continue
            if cr > 0 and left is None:
                left = d
            elif cr < 0 and right is None:
                right = d
        if left is None or right is None or abs(left - right) > 0.3:
            return None
        return c, left, right, S

    @staticmethod
    def _wall_ahead(S, c: Pt, nt: Pt, r: float) -> bool:
        """Круг радиуса r с центром c упёрся в контур впереди (развилка или торец)."""
        for d, q, _e in S:
            if d > 1.02 * r + 0.2:
                break
            if vdot(vsub(q, c), nt) > 0.75 * d:
                return True
        return False

    def trace_medial(self, start: Pt, t: Pt, w: float, maxlen: float, stop=None):
        """Ведёт ось от start по направлению t по настоящей середине контура.
        Возвращает (точки, причина, радиусы): 'fork' (впереди стенка: стык или торец), 'stop'
        (сработал stop(c, касательная)), 'len' (длина исчерпана) или 'fail' (середина потеряна).
        Радиус вписанного круга в каждой точке даёт настоящую толщину штриха в стыке."""
        cap = 1.7 * w
        h = max(1.5, 0.05 * w)
        s0 = self._snap(start, t, cap)
        if s0 is None:
            return [], "fail", []
        c, dl, dr, _S = s0
        r = 0.5 * (dl + dr)
        pts: List[Pt] = [c]
        rad: List[float] = [r]
        length = 0.0
        while length < maxlen:
            pred = (c[0] + t[0] * h, c[1] + t[1] * h)
            s1 = self._snap(pred, t, cap)
            if s1 is None:
                return pts, "fail", rad
            c1, dl, dr, S = s1
            nt = vunit(vsub(c1, c))
            r1 = 0.5 * (dl + dr)
            if nt is None or vdot(nt, t) < 0.7 or abs(r1 - r) > 1.3 * h or not self.inside(c1):
                return pts, "fail", rad  # перескочили на другую ветвь или вышли из материала
            if length >= 0.6 * r and self._wall_ahead(S, c1, nt, r1):
                lo, hi, r_hi = c, c1, r1  # уточняем место, где круг касается стенки
                for _ in range(6):
                    sm = self._snap(lerp(lo, hi, 0.5), t, cap)
                    if sm is None:
                        break
                    if self._wall_ahead(sm[3], sm[0], nt, 0.5 * (sm[1] + sm[2])):
                        hi, r_hi = sm[0], 0.5 * (sm[1] + sm[2])
                    else:
                        lo = sm[0]
                pts.append(hi)
                rad.append(r_hi)
                return pts, "fork", rad
            pts.append(c1)
            rad.append(r1)
            length += dist(c, c1)
            c, t, r = c1, nt, r1
            if stop is not None and stop(c, t):
                return pts, "stop", rad
        return pts, "len", rad


    # --- продолжение кривой чаши до оси стойки ---------------------------------------
    def _inward_points(self, e: "End", W: float, ends: Optional[List["End"]]) -> List[Pt]:
        """Точки оси от конца e вглубь штриха на длину W. Короткий штрих (кусок чаши между узлами контура)
        продолжается на соседний, если стык с ним почти гладкий."""
        pts: List[Pt] = [e.pos]
        length = 0.0
        cur = e
        seen = {id(e.s)}
        for _hop in range(4):
            poly = self.polys[id(cur.s)]
            seq = list(reversed(poly)) if cur.idx == 1 else list(poly)
            for q in seq[1:]:
                length += dist(pts[-1], q)
                pts.append(q)
                if length >= W:
                    return pts
            far = next((x for x in (ends or []) if x.s is cur.s and x.idx != cur.idx), None)
            if far is None:
                break
            nxt, best = None, 1e18
            for x in ends or []:
                if x.s is cur.s or id(x.s) in seen or x.used:
                    continue
                g = dist(x.pos, far.pos)
                if g > 0.7 * max(far.w, x.w) or abs(angle_deg(far.t, (-x.t[0], -x.t[1]))) > 35.0:
                    continue
                if g > 1e-6 and vdot(far.t, vsub(x.pos, far.pos)) < -0.1 * g:
                    continue
                if g < best:
                    nxt, best = x, g
            if nxt is None:
                break
            seen.add(id(nxt.s))
            cur = nxt
        return pts

    def continue_arm(self, e: "End", s2: "Stroke", ends: Optional[List["End"]]) -> Optional[_Link]:
        key = ("cont", id(e), id(s2))
        if key not in self._memo:
            self._memo[key] = self._continue_arm(e, s2, ends)
        return self._memo[key]

    def _continue_arm(self, e: "End", s2: "Stroke", ends: Optional[List["End"]]) -> Optional[_Link]:
        """Чаша b, d, p, q: от вершины чаши (где ось идёт перпендикулярно стойке) ось продолжается одной
        дугой окружности и садится на ось стойки под углом prm.arm_angle. Участок оси от вершины до
        конца дуги заменяется: у выреза измеренная ось закручивается внутрь (вырез тянет середину
        между контурами), а чаша, как её рисует человек, остаётся овалом. Замена допускается,
        пока дуга отходит от измеренной оси не больше чем на 0.2 толщины."""
        prm = self.prm
        poly = self.polys[id(s2)]
        if dist_to_polyline(e.pos, poly) < 0.25 * e.w:
            return None
        A1 = poly[0]
        ul = vunit(vsub(poly[-1], poly[0]))
        if ul is None:
            return None
        nl = (-ul[1], ul[0])
        Ls = dist(poly[0], poly[-1])
        own = self.polys[id(e.s)]
        chain = self._inward_points(e, 4.0 * e.w, ends)
        if len(chain) < 8:
            return None
        # вершина чаши: первая точка пути от конца, где ход оси становится перпендикулярным стойке
        g = [vdot(vsub(chain[i + 1], chain[i]), ul) for i in range(len(chain) - 1)]
        apex = None
        for i in range(1, len(g)):
            if g[i] * g[0] <= 0.0:
                apex = i
                break
        if apex is None:
            return None
        Q = chain[apex]
        if dist_to_polyline(Q, own) > 0.5:
            return None  # вершина на соседнем штрихе (чаша разорвана): отрезать нечего
        old = chain[:apex + 1]  # заменяемый участок: от конца дуги до вершины
        if sum(dist(a_, b_) for a_, b_ in zip(old, old[1:])) > 2.5 * e.w:
            return None
        tq = vunit(vsub(chain[apex - 1], chain[min(apex + 1, len(chain) - 1)]))  # ход в сторону конца дуги
        if tq is None:
            return None
        # к центру кривизны: сторона, куда изгибается чаша за вершиной
        n0 = (-tq[1], tq[0])
        far = chain[min(len(chain) - 1, apex + 6)]
        n = n0 if vdot(vsub(far, Q), n0) >= 0 else (-n0[0], -n0[1])
        sigma = 1 if vcross(tq, n) > 0 else -1
        sdq = vdot(vsub(Q, A1), nl)
        # угол поворота: на оси стойки ход составляет с ней угол arm_angle
        target = math.sin(math.radians(max(10.0, min(80.0, prm.arm_angle))))

        def f_inc(th: float) -> float:
            h = (tq[0] * math.cos(th) + n[0] * math.sin(th), tq[1] * math.cos(th) + n[1] * math.sin(th))
            return abs(vdot(h, nl)) - target

        lo, hi = 0.0, math.radians(100.0)
        if f_inc(lo) <= 0.0 or f_inc(hi) >= 0.0:
            return None  # вершина уже положе нужного угла или повернуть до него нельзя
        for _ in range(40):
            mid = 0.5 * (lo + hi)
            if f_inc(mid) > 0.0:
                lo = mid
            else:
                hi = mid
        theta = 0.5 * (lo + hi)
        den = math.sin(theta) * vdot(tq, nl) + (1.0 - math.cos(theta)) * vdot(n, nl)
        if abs(den) < 1e-9 or -sdq / den <= 0.0:
            return None
        rho = -sdq / den
        if not (0.5 * e.w <= rho <= 25.0 * e.w):
            return None
        C = (Q[0] + n[0] * rho, Q[1] + n[1] * rho)
        cubics = _arc_cubics(C, rho, Q, tq, sigma, theta)
        X = (Q[0] + rho * (math.sin(theta) * tq[0] + (1 - math.cos(theta)) * n[0]),
             Q[1] + rho * (math.sin(theta) * tq[1] + (1 - math.cos(theta)) * n[1]))
        cubics[-1] = (cubics[-1][0], cubics[-1][1], cubics[-1][2], X)
        al = vdot(vsub(X, A1), ul)
        if not (-0.02 * Ls <= al <= 1.02 * Ls):
            return None
        w2 = _width_at(s2, X)
        if rho * theta > 3.0 * max(e.w, w2):
            return None
        samples = [q for c in cubics for q in flatten(c, 16)]
        # дуга не должна уходить от измеренной оси (на заменяемом участке) дальше 0.15 толщины
        if max(dist_to_polyline(q, samples) for q in old) > 0.2 * e.w:
            return None
        if not all(self.inside(q) for q in samples[1:-1]):
            return None
        return _Link(cubics, "mid", [(q, e.w) for q in samples[1::2]], cut=Q, how="oval")

    def medial_arm(self, e: "End", s2: "Stroke") -> Optional[_Link]:
        """Конец дуги e доводится до оси прямой стойки s2 по середине контура. Ось дуги идёт
        вдоль настоящей срединной линии между вырезом и чашей, затем по касательной входит в
        ось стойки. Если середина теряется, возвращает None (останутся прежние способы)."""
        key = ("arm", id(e), id(s2))
        if key not in self._memo:
            self._memo[key] = self._medial_arm(e, s2)
        return self._memo[key]

    def _medial_arm(self, e: "End", s2: "Stroke") -> Optional[_Link]:
        poly = self.polys[id(s2)]
        if dist_to_polyline(e.pos, poly) < 0.25 * e.w:
            return None  # конец уже лежит на оси цели (дуга продолжает стойку): это G1, а не T
        w2 = _width_at(s2, e.pos)
        wmax = max(e.w, w2)
        pts, why, rad = self.trace_medial(e.pos, e.t, wmax, 3.0 * wmax,
                                          lambda c, nt: dist_to_polyline(c, poly) <= 0.25)
        if why not in ("fork", "stop") or len(pts) < 2:
            return None
        if dist(pts[0], e.pos) > 0.2 * e.w:
            return None  # ось дуги у конца ушла от середины: доверять трассе нельзя
        F = pts[-1]
        dF = dist_to_polyline(F, poly)
        if dF > 0.4 * w2 + 1.0:
            return None  # круг упёрся в стенку не у этой стойки
        X = F
        if dF > 0.25:
            tF = vunit(vsub(pts[-1], pts[max(0, len(pts) - 3)])) or e.t
            best = None
            for a_, b_ in zip(poly, poly[1:]):  # продолжение по касательной до оси стойки
                r_ = _seg_ray(F, tF, a_, b_)
                if r_ is not None and 0.0 <= r_[1] <= 1.0 and r_[0] >= 0.0 and (best is None or r_[0] < best[0]):
                    best = r_
            if best is not None and best[0] <= 0.5 * w2 + 1.0:
                X = (F[0] + tF[0] * best[0], F[1] + tF[1] * best[0])
            else:
                X = min((_nearest_on_segment(F, a_, b_) for a_, b_ in zip(poly, poly[1:])),
                        key=lambda q: dist(q, F))
            if dist(F, X) > 0.5 * w2 or not self.seg_inside(F, X, 3):
                return None
        path = [e.pos] + pts[1:] + ([X] if dist(X, F) > 0.05 else [])
        cubics = _fit_path(path, e.t, None, max(0.35, 0.0005 * self.upm))
        if not cubics or not all(self.inside(q) for c in cubics for q in flatten(c, 10)[1:-1]):
            return None
        prof = [(q, 2.0 * r) for q, r in zip(pts[1:], rad[1:])]
        if dist(X, F) > 0.05:
            prof.append((X, 2.0 * rad[-1]))
        return _Link(cubics, "mid", prof, how="medial")

    def _fork_partner(self, e1: "End", e2: "End", X: Pt, ends: Optional[List["End"]]) -> bool:
        """Есть ли у узла L третий свободный конец (перекладина e, ножка k): тогда узел остаётся
        развилкой Y, а не гладким мостом между двумя дугами."""
        for e3 in ends or []:
            if e3.used or e3.terminal or e3.s in (e1.s, e2.s):
                continue
            g = vsub(X, e3.pos)
            along, perp = vdot(g, e3.t), abs(vcross(g, e3.t))
            if -0.3 * e3.w <= along <= self.prm.ext_factor * max(e3.w, e1.w, e2.w) / 0.5 and perp <= 0.4 * e3.w:
                return True
        return False

    def medial_bridge(self, e1: "End", e2: "End") -> Optional[_Link]:
        """Два конца дуг одной чаши (разрыв у узла контура): ось ведётся от e1 по середине контура
        до e2. Это замена угла L и «усов» G1, когда обе дуги лежат на одной срединной линии."""
        key = ("bridge", id(e1), id(e2))
        if key not in self._memo:
            self._memo[key] = self._medial_bridge(e1, e2)
        return self._memo[key]

    def _medial_bridge(self, e1: "End", e2: "End") -> Optional[_Link]:
        wmax = max(e1.w, e2.w)
        target = e2.pos
        h = max(1.5, 0.05 * wmax)

        def stop(c, nt):
            return dist(c, target) <= 1.5 * h or vdot(vsub(target, c), nt) <= 0.0

        pts, why, rad = self.trace_medial(e1.pos, e1.t, wmax, 2.5 * wmax + dist(e1.pos, e2.pos), stop)
        if why != "stop" or len(pts) < 2:
            return None
        if dist(pts[0], e1.pos) > 0.2 * e1.w or dist(pts[-1], target) > 0.35 * wmax:
            return None
        tin = vunit(vsub(target, pts[-1])) or e2.t
        if abs(angle_deg(tin, (-e2.t[0], -e2.t[1]))) > 30.0:
            return None  # подход к e2 не вдоль её оси: это другая ветвь
        cubics = _fit_path([e1.pos] + pts[1:] + [target], e1.t, (-e2.t[0], -e2.t[1]),
                           max(0.35, 0.0005 * self.upm))
        if not cubics or not all(self.inside(q) for c in cubics for q in flatten(c, 10)[1:-1]):
            return None
        return _Link(cubics, "mid", [(q, 2.0 * r) for q, r in zip(pts[1:], rad[1:])])

    @staticmethod
    def trim_axis(e: "End", Q: Pt) -> None:
        """Обрезает ось штриха с конца e в точке Q (деление кубики де Кастельжо)."""
        axis, tags = e.s.axis, e.s.tags
        best = None
        for i, c in enumerate(axis):
            for k in range(0, 49):
                u = k / 48.0
                d = dist(bez_pt(c, u), Q)
                if best is None or d < best[0]:
                    best = (d, i, u)
        d, i, u = best
        c = axis[i]
        # де Кастельжо: левая (0..u) и правая (u..1) части
        p01, p12, p23 = lerp(c[0], c[1], u), lerp(c[1], c[2], u), lerp(c[2], c[3], u)
        p012, p123 = lerp(p01, p12, u), lerp(p12, p23, u)
        m_ = lerp(p012, p123, u)
        left, right = (c[0], p01, p012, m_), (m_, p123, p23, c[3])
        def tiny(q: Cubic) -> bool:
            return dist(q[0], q[3]) < 0.5 and dist(q[0], q[1]) < 0.5 and dist(q[0], q[2]) < 0.5

        if e.idx == 1:
            keep = [] if tiny(left) else [left]
            e.s.axis = axis[:i] + keep
            e.s.tags = tags[:i] + [tags[i]] * len(keep)
        else:
            keep = [] if tiny(right) else [right]
            e.s.axis = keep + axis[i + 1:]
            e.s.tags = [tags[i]] * len(keep) + tags[i + 1:]

    @staticmethod
    def apply_link(e: "End", link: _Link) -> None:
        """Приставляет участок link (от e.pos наружу) к оси штриха на нужном конце.
        Толщина вдоль участка (радиус вписанного круга) дописывается в профиль штриха."""
        e.s.prof.extend(link.prof)
        if link.cut is not None:
            _Connector.trim_axis(e, link.cut)
        if e.idx == 1:
            for c in link.cubics:
                e.s.axis.append(c)
                e.s.tags.append(link.tag)
        else:
            for c in link.cubics:
                e.s.axis.insert(0, _reverse_cubic(c))
                e.s.tags.insert(0, link.tag)

    def remaining(self, e: End) -> float:
        """Сколько самой длинной стороны штриха осталось за концом оси. Если ребро-сторона идёт
        дальше конца оси (стойка b, d, p, q ниже чаши), штрих идёт дальше, и угол L с ним
        неверен: это T на сквозной штрих."""
        best = None
        for p in e.s.pieces:
            if p.kind == "blob":
                continue
            for k, pt in enumerate((p.axis[0], p.axis[3])):
                d = dist(pt, e.pos)
                if best is None or d < best[0]:
                    best = (d, p, k)
        if best is None:
            return 0.0
        _, p, k = best
        sa = p.samples
        ta0, ta1, tb0, tb1 = sa[0][0], sa[-1][0], sa[0][1], sa[-1][1]
        fwd = tb1 >= tb0  # B идёт в ту же сторону по t, что и A
        if k == 1:
            ra = (1.0 - ta1) * p.a.length
            rb = (1.0 - tb1) * p.b.length if fwd else tb1 * p.b.length
        else:
            ra = ta0 * p.a.length
            rb = tb0 * p.b.length if fwd else (1.0 - tb0) * p.b.length
        return max(ra, rb)

    def end_curvature(self, e: End) -> float:
        """Знаковая кривизна оси у конца e в направлении наружу (плюс: поворот против часовой)."""
        idxs = range(len(e.s.axis) - 1, -1, -1) if e.idx == 1 else range(len(e.s.axis))
        for i in idxs:
            if e.s.tags[i] != "fit":
                continue
            P0, P1, P2, P3 = e.s.axis[i]
            if e.idx == 1:
                d1 = (3 * (P3[0] - P2[0]), 3 * (P3[1] - P2[1]))
                d2 = (6 * (P3[0] - 2 * P2[0] + P1[0]), 6 * (P3[1] - 2 * P2[1] + P1[1]))
                sign = 1.0
            else:
                d1 = (3 * (P1[0] - P0[0]), 3 * (P1[1] - P0[1]))
                d2 = (6 * (P2[0] - 2 * P1[0] + P0[0]), 6 * (P2[1] - 2 * P1[1] + P0[1]))
                sign = -1.0  # наружу значит назад по кубике
            n1 = math.hypot(d1[0], d1[1])
            return sign * vcross(d1, d2) / n1 ** 3 if n1 > 1e-9 else 0.0
        return 0.0

    def tangent_merge(self, e: End, target: "Stroke") -> Optional[Cubic]:
        """Дуга входит в прямую стойку по касательной: кривая идёт от конца дуги к основанию
        перпендикуляра на ось стойки и приходит туда вдоль стойки, в сторону поворота дуги.
        Так чаша b, d, p, q покидает стойку вверху и внизу, а стойка между ними служит левой
        стороной кольца."""
        poly = self.polys[id(target)]
        best = None
        for a_, b_ in zip(poly, poly[1:]):
            q = _nearest_on_segment(e.pos, a_, b_)
            dd = dist(q, e.pos)
            if best is None or dd < best[0]:
                best = (dd, q, vunit(vsub(b_, a_)))
        if best is None or best[2] is None:
            return None
        D, F, u = best
        w2 = _width_at(target, F)
        if D > 0.5 * w2 + 1.5 * self.prm.ext_factor * e.w or abs(vcross(e.t, u)) < 0.2:
            return None
        k = self.end_curvature(e)
        if abs(k) * e.w > 1e-3:
            d_t = u if vcross(e.t, u) * k > 0 else (-u[0], -u[1])
        else:
            d_t = u if vdot(u, e.t) >= 0 else (-u[0], -u[1])
        # настоящая дуга окружности: поворот на угол theta до параллели со стойкой на расстоянии D,
        # радиус R = D / (1 - cos theta); точка входа смещена вдоль стойки на R sin theta,
        # иначе кривой пришлось бы сначала спуститься, а потом вернуться (волна)
        cos_t = max(-1.0, min(1.0, vdot(e.t, d_t)))
        theta = math.acos(cos_t)
        # Круговая конструкция устойчива только в умеренном диапазоне углов.
        # Вне него вызывающий код пробует corridor_path(), а не прямой fallback.
        if theta < math.radians(8) or theta > math.radians(40):
            return None
        R = D / (1.0 - cos_t)
        sh = R * math.sin(theta)
        exact = sh <= 1.5 * D + 0.3 * e.w
        sh = min(sh, 1.5 * D + 0.3 * e.w)  # конец, идущий почти вдоль стойки, сходится к ней плавно, а не по идеальной окружности
        F2 = (F[0] + d_t[0] * sh, F[1] + d_t[1] * sh)
        if dist_to_polyline(F2, poly) > 0.15 * w2:
            return None  # точка входа вышла за ось стойки
        h = 4.0 / 3.0 * math.tan(theta / 4.0) * R if exact else 0.4 * dist(e.pos, F2)
        c: Cubic = (e.pos, (e.pos[0] + e.t[0] * h, e.pos[1] + e.t[1] * h),
                    (F2[0] - d_t[0] * h, F2[1] - d_t[1] * h), F2)
        if not all(self.inside(q) for q in flatten(c, 10)[1:-1]):
            return None
        return c

    def foot_link(self, e: End, target: "Stroke") -> Optional[Cubic]:
        """Кратчайшая прямая от конца до оси цели (основание перпендикуляра), если она лежит в материале."""
        poly = self.polys[id(target)]
        best = None
        for a_, b_ in zip(poly, poly[1:]):
            q = _nearest_on_segment(e.pos, a_, b_)
            dd = dist(q, e.pos)
            if best is None or dd < best[0]:
                best = (dd, q)
        if best is None:
            return None
        D, F = best
        if D < 0.5 or D > 0.5 * _width_at(target, F) + 1.5 * self.prm.ext_factor * e.w:
            return None
        if not self.seg_inside(e.pos, F, 4):
            return None
        return line_cubic(e.pos, F)

    def through_extend(self, e: End, rem: float = 1e9) -> bool:
        """Сквозной штрих: ось достраивается до торца, который виден в контуре, но не дальше, чем
        продолжается его сторона (остаток ребра rem) плюс полтолщины."""
        h = self.ahead(e.pos, e.t)
        if h is None or h[0] < 0.05 * e.w or h[0] > max(rem + 0.5 * e.w, 4.0 * e.w):
            return False
        n = h[1].normal(h[2])
        if -(e.t[0] * n[0] + e.t[1] * n[1]) < math.cos(math.radians(50)) or not self.seg_inside(e.pos, h[3], 4):
            return False
        self.extend(e, h[3], "cap")
        self.record_cap(e, h)
        e.pos = h[3]
        e.terminal = True
        e.mark = "terminal"
        return True

    # --- шаги -----------------------------------------------------------------------
    def bridge_gaps(self):
        """Щели внутри цепочки осей (стык с другим штрихом) закрываются плавным мостом."""
        for s in self.strokes:
            n = len(s.axis)
            if n < 2 and not s.closed:
                continue
            axis: List[Cubic] = []
            tags: List[str] = []
            for j in range(n):
                axis.append(s.axis[j])
                tags.append(s.tags[j])
                if j == n - 1 and not s.closed:
                    break
                nxt = s.axis[(j + 1) % n]
                gap = dist(s.axis[j][3], nxt[0])
                chord = vsub(nxt[0], s.axis[j][3])
                ta = end_tangents(s.axis[j])[1]
                tb = end_tangents(nxt)[0]
                # мост короче 1.5 толщины, а его «усы» смотрят вперёд: иначе ось заворачивается назад
                if 0.5 < gap <= 1.5 * s.w_mean and vdot(ta, chord) > 0.1 * gap and vdot(tb, chord) > 0.1 * gap \
                        and self.seg_inside(s.axis[j][3], nxt[0]):
                    h = gap / 3.0
                    axis.append((s.axis[j][3],
                                 (s.axis[j][3][0] + ta[0] * h, s.axis[j][3][1] + ta[1] * h),
                                 (nxt[0][0] - tb[0] * h, nxt[0][1] - tb[1] * h), nxt[0]))
                    tags.append("bridge")
            s.axis, s.tags = axis, tags

    def make_ends(self) -> List[End]:
        ends: List[End] = []
        for s in self.strokes:
            if s.closed:
                continue
            c0, c1 = s.axis[0], s.axis[-1]
            t0 = end_tangents(c0)[0]
            t1 = end_tangents(c1)[1]
            ends.append(End(s, 0, c0[0], (-t0[0], -t0[1]), _width_at(s, c0[0])))
            ends.append(End(s, 1, c1[3], t1, _width_at(s, c1[3])))
        for e in ends:
            h = self.ahead(e.pos, e.t)
            e.terminal = h is not None and h[0] <= self.prm.term_factor * e.w
            if e.terminal:
                self.record_cap(e, h)
        return ends

    def extend(self, e: End, X: Pt, tag: str = "ext"):
        if dist(e.pos, X) < 0.5:
            return
        if e.idx == 1:
            e.s.axis.append(line_cubic(e.pos, X))
            e.s.tags.append(tag)
        else:
            e.s.axis.insert(0, line_cubic(X, e.pos))
            e.s.tags.insert(0, tag)

    def candidates(self, ends: List[End]) -> list:
        """Кандидаты в соединения: (вид, цена, конец1, конец2, точка, штрих-цель, толщина, sin угла)."""
        prm = self.prm
        cands = []
        free = [e for e in ends if not e.terminal]
        for e in free:  # продление до оси другого штриха (T)
            for s2 in self.strokes:
                if s2 is e.s:
                    continue
                best = None
                poly = self.polys[id(s2)]
                for a_, b_ in zip(poly, poly[1:]):
                    r = _seg_ray(e.pos, e.t, a_, b_)
                    if r is None or not (0.0 <= r[1] <= 1.0) or r[0] < -0.3 * e.w:
                        continue
                    u = vunit(vsub(b_, a_))
                    sin_t = abs(vcross(e.t, u)) if u else 1.0
                    X = (e.pos[0] + e.t[0] * r[0], e.pos[1] + e.t[1] * r[0])
                    w2 = _width_at(s2, X)
                    # конец штриха упирается в край цели (полтолщины) плюс непокрытый участок стыка
                    if r[0] <= (0.5 * w2 + prm.ext_factor * e.w) / max(sin_t, 0.5) and (best is None or r[0] < best[0]):
                        best = (r[0], X, w2, sin_t)
                if best:
                    sd, X, w2, sin_t = best
                    cands.append(("T", max(sd, 0.0) / max(e.w, w2), e, None, X, s2, w2, sin_t))
        for i in range(len(free)):  # два конца: угол (L) или плавное продолжение (G1)
            for j in range(i + 1, len(free)):
                e1, e2 = free[i], free[j]
                if e1.s is e2.s:
                    continue
                wm = max(e1.w, e2.w)
                g = vsub(e2.pos, e1.pos)
                ang = abs(angle_deg(e1.t, (-e2.t[0], -e2.t[1])))
                if ang <= self.G1_DEG:
                    along, perp = vdot(g, e1.t), abs(vcross(g, e1.t))
                    if -0.3 * wm <= along <= prm.ext_factor * wm and perp <= 0.4 * wm:
                        cands.append(("G1", (abs(along) + perp) / wm, e1, e2, lerp(e1.pos, e2.pos, 0.5), None, wm, 1.0))
                    continue
                cr = vcross(e1.t, e2.t)
                if abs(cr) < 1e-6:
                    continue
                s1, s2 = vcross(g, e2.t) / cr, vcross(g, e1.t) / cr
                k = max(abs(cr), 0.35)
                if s1 < -0.3 * e1.w or s2 < -0.3 * e2.w:
                    continue
                if s1 > (0.5 * e2.w + prm.ext_factor * e1.w) / k or s2 > (0.5 * e1.w + prm.ext_factor * e2.w) / k:
                    continue
                X = (e1.pos[0] + e1.t[0] * s1, e1.pos[1] + e1.t[1] * s1)
                cands.append(("L", max(s1, 0) / e1.w + max(s2, 0) / e2.w, e1, e2, X, None, wm, abs(cr)))
        cands.sort(key=lambda c: c[1])
        return cands

    def plan(self, cands: list, banned: set, ends: Optional[List[End]] = None) -> list:
        """Выбор соединений без изменения осей. Дешёвые идут первыми; конец служит одному соединению."""
        prm = self.prm
        used: set = set()
        res = []
        for c in cands:
            kind, cost, e1, e2, X, s2, wm, sin_t = c
            if id(e1) in used or (e2 is not None and id(e2) in used) or (id(e1), id(e2)) in banned:
                continue
            sup = self.support(X, e1.w, wm if e2 is None else e2.w, sin_t)
            path: Optional[_Link] = None
            if kind in ("L", "G1") and prm.medial_joints and e2 is not None \
                    and e1.s.kind == "arc" and e2.s.kind == "arc" \
                    and not self._fork_partner(e1, e2, X, ends):
                # две дуги одной чаши, разорванные у узла контура: ось по середине контура вместо «уса»
                path = self.medial_bridge(e1, e2)
                if path is not None:
                    kind, X = "G1", lerp(e1.pos, e2.pos, 0.5)
            if kind != "G1" and sup < prm.min_support:
                continue  # нет опоры в существующем контуре (у гладкого продолжения G1 углов нет)
            if kind == "T":
                # Т-стык с участием дуги должен быть плавным. Раньше прямая e1->X
                # проходила проверку seg_inside(), после чего extend() добавлял её
                # даже при неудаче tangent_merge(). Отсюда появлялись прямые участки
                # между чашей и стойкой в b/d/p/q.
                curved_join = e1.s.kind == "arc" or s2.kind == "arc"
                arc_to_stem = e1.s.kind == "arc" and s2.kind == "line"
                if arc_to_stem and prm.arm_continue:
                    path = self.continue_arm(e1, s2, ends)  # продолжение кривой чаши
                if path is None and arc_to_stem and prm.medial_joints:
                    path = self.medial_arm(e1, s2)  # ось по настоящей середине контура
                if path is None and arc_to_stem:
                    cub = self.tangent_merge(e1, s2)
                    path = _Link([cub]) if cub is not None else None
                if path is None and curved_join:
                    cub = self.corridor_path(e1, X, s2)
                    path = _Link([cub]) if cub is not None else None
                if curved_join:
                    # Не подменяем неудачную кривую прямым отрезком. Лучше оставить
                    # конец открытым, чем экспортировать геометрически неверное ребро.
                    if path is None:
                        continue
                    X = path.end
                elif not self.seg_inside(e1.pos, X, 3):
                    # Для прямого штриха прямое соединение допустимо только если
                    # оно лежит в материале; иначе ищем центральную линию прохода.
                    cub = self.foot_link(e1, s2) or self.corridor_path(e1, X, s2)
                    if cub is None:
                        continue
                    path = _Link([cub])
                    X = path.end
            elif path is not None:
                pass
            elif kind == "G1":
                if not self.seg_inside(e1.pos, e2.pos):
                    continue
            elif not self.seg_inside(e1.pos, X, 3) or (e2 is not None and not self.seg_inside(e2.pos, X, 3)):
                continue
            if kind != "G1" and path is None and not self.near_ink(X, 0.15 * min(e1.w, wm)):
                continue  # точка стыка вне буквы
            used.add(id(e1))
            if e2 is not None:
                used.add(id(e2))
            res.append((kind, cost, e1, e2, X, s2, wm, sin_t, sup, path))
        # два рукава на одну стойку (верх и низ чаши) строятся одним способом: иначе чаша несимметрична
        by_stem: Dict[int, List[int]] = {}
        for i, r in enumerate(res):
            if r[0] == "T" and r[9] is not None and r[9].how in ("oval", "medial"):
                by_stem.setdefault(id(r[5]), []).append(i)
        for idxs in by_stem.values():
            if len({res[i][9].how for i in idxs}) > 1:
                for i in idxs:
                    kind, cost, e1, e2, X, s2, wm, sin_t, sup, path = res[i]
                    if path.how == "oval":
                        alt = self.medial_arm(e1, s2)
                        if alt is not None:
                            res[i] = (kind, cost, e1, e2, alt.end, s2, wm, sin_t, sup, alt)
        return res

    def reroute(self, plan: list, cands: list) -> set:
        """Угол L двух штрихов, рядом с которым проходит ось третьего штриха той же части
        буквы (рука и нога k рядом со стойкой), заменяется двумя T на этот третий штрих."""
        banned = set()
        has_T = {(id(c[2]), id(c[5])) for c in cands if c[0] == "T"}
        for kind, cost, e1, e2, X, s2, wm, sin_t, sup, path in plan:
            if kind != "L":
                continue
            for s3 in self.strokes:
                if s3 is e1.s or s3 is e2.s or not self.same_ink(e1.s, s3):
                    continue
                d3 = dist_to_polyline(X, self.polys[id(s3)])
                w3 = _width_at(s3, X)
                if d3 > 0.5 * w3 and (id(e1), id(s3)) in has_T and (id(e2), id(s3)) in has_T:
                    banned.add((id(e1), id(e2)))
                    break
        return banned

    def run(self):
        prm = self.prm
        self.bridge_gaps()
        self.polys = {id(s): _axis_poly(s) for s in self.strokes}
        ends = self.make_ends()
        for _ in range(2):  # сквозные штрихи: угол L с ними заменяется на T
            changed = False
            for c in self.candidates(ends):
                if c[0] != "L":
                    continue
                e1, e2 = c[2], c[3]
                r1, r2 = self.remaining(e1), self.remaining(e2)
                t1 = r1 > 1.0 * e2.w + 0.35 * e1.w
                t2 = r2 > 1.0 * e1.w + 0.35 * e2.w
                for e, mine, other, rem in ((e1, t1, t2, r1), (e2, t2, t1, r2)):
                    if mine and not other and not e.terminal and self.through_extend(e, rem):
                        changed = True
            if not changed:
                break
            self.polys = {id(s): _axis_poly(s) for s in self.strokes}
        cands = self.candidates(ends)
        banned: set = set()
        plan = self.plan(cands, banned, ends)
        for _ in range(3):
            more = self.reroute(plan, cands) - banned
            if not more:
                break
            banned |= more
            plan = self.plan(cands, banned, ends)
        for kind, cost, e1, e2, X, s2, wm, sin_t, sup, path in plan:
            if kind == "G1":
                if path is not None:  # ось по середине контура между двумя дугами
                    self.apply_link(e1, path)
                else:
                    h = dist(e1.pos, e2.pos) / 3.0
                    c = (e1.pos, (e1.pos[0] + e1.t[0] * h, e1.pos[1] + e1.t[1] * h),
                         (e2.pos[0] + e2.t[0] * h, e2.pos[1] + e2.t[1] * h), e2.pos)
                    self.apply_link(e1, _Link([c], "bridge"))
                e1.used = e2.used = True
                self.joints.append(Joint("G1", X, [e1.s, e2.s], sup))
                continue
            if path is not None:
                self.apply_link(e1, path)
            else:
                self.extend(e1, X)
            e1.used = True
            if kind == "L":
                self.extend(e2, X)
                e2.used = True
                j = Joint("L", X, [e1.s, e2.s], sup)
                # третий свободный конец, чья линия проходит через этот же узел (развилка Y)
                for e3 in ends:
                    if e3.used or e3.terminal or e3.s in (e1.s, e2.s):
                        continue
                    g = vsub(X, e3.pos)
                    along, perp = vdot(g, e3.t), abs(vcross(g, e3.t))
                    if -0.3 * e3.w <= along <= prm.ext_factor * max(e3.w, e1.w, e2.w) / 0.5 and perp <= 0.4 * e3.w \
                            and self.seg_inside(e3.pos, X, 3):
                        self.extend(e3, X)
                        e3.used = True
                        j.strokes.append(e3.s)
                        j.kind = "Y"
                if j.kind == "L":
                    self.corners.append((X, e1, e2))
                self.joints.append(j)
            else:
                self.joints.append(Joint("T", X, [e1.s, s2], sup))
        # свободные концы: достраиваем до торца, который виден в контуре
        cos_cap = math.cos(math.radians(50))
        for e in ends:
            if e.used:
                continue
            if e.terminal:
                e.mark = "terminal"
                continue
            h = self.ahead(e.pos, e.t)
            if h is not None and prm.cap_factor > 0 and 0.05 * e.w < h[0] <= prm.cap_factor * e.w:
                n = h[1].normal(h[2])
                near_convex = any(v.kind == "convex" and dist(v.pos, h[3]) <= 0.8 * e.w for v in self.a.vertices)
                if -(e.t[0] * n[0] + e.t[1] * n[1]) >= cos_cap and near_convex:
                    self.extend(e, h[3], "cap")
                    self.record_cap(e, h)
                    e.mark = "cap"
                    continue
            e.mark = "open"
        self.polys = {id(s): _axis_poly(s) for s in self.strokes}
        self.crossings()
        for s in self.strokes:
            s.length = sum(_cubic_len(c) for c in s.axis)
        self.a.joints = self.joints
        self.a.end_marks = [(e.pos, e.mark) for e in ends if e.mark in ("open", "cap")]

    def crossings(self):
        """Оси двух штрихов пересекаются вдали от концов: узел X (скажем, диагонали x)."""
        for i, s1 in enumerate(self.strokes):
            p1 = self.polys[id(s1)]
            b1 = (min(q[0] for q in p1), min(q[1] for q in p1), max(q[0] for q in p1), max(q[1] for q in p1))
            for s2 in self.strokes[i + 1:]:
                p2 = self.polys[id(s2)]
                b2 = (min(q[0] for q in p2), min(q[1] for q in p2), max(q[0] for q in p2), max(q[1] for q in p2))
                if b1[0] > b2[2] or b2[0] > b1[2] or b1[1] > b2[3] or b2[1] > b1[3]:
                    continue
                for a_, b_ in zip(p1, p1[1:]):
                    d = vsub(b_, a_)
                    L = math.hypot(d[0], d[1])
                    if L < 1e-9:
                        continue
                    u1 = (d[0] / L, d[1] / L)
                    for c_, e_ in zip(p2, p2[1:]):
                        r = _seg_ray(a_, u1, c_, e_)
                        if r is None or not (0.0 <= r[0] <= L) or not (0.0 <= r[1] <= 1.0):
                            continue
                        X = (a_[0] + u1[0] * r[0], a_[1] + u1[1] * r[0])
                        w1, w2 = _width_at(s1, X), _width_at(s2, X)
                        lim = 0.5 * max(w1, w2)
                        if any(dist(X, q) < lim for q in (p1[0], p1[-1], p2[0], p2[-1])):
                            continue
                        if any(dist(X, j.pos) < lim for j in self.joints):
                            continue
                        self.joints.append(Joint("X", X, [s1, s2], self.support(X, w1, w2)))


def _outside_ink(a: "Analysis", pt: Pt, tol: float) -> bool:
    """Точка вне материала и дальше tol от контура."""
    if winding_number(pt, a.polys) != 0:
        return False
    for e in a.edges:
        bb = e.bbox
        if bb[0] - tol <= pt[0] <= bb[2] + tol and bb[1] - tol <= pt[1] <= bb[3] + tol \
                and dist_to_polyline(pt, e.poly) <= tol:
            return False
    return True


def _clamp_to_ink(a: "Analysis", side: List[Pt], end_i: int, prev_i: int) -> None:
    """Угол ленты, выступивший за контур, возвращается по своей стороне к границе материала."""
    c, p0 = side[end_i], side[prev_i]
    if not _outside_ink(a, c, 1.0) or winding_number(p0, a.polys) == 0:
        return
    lo, hi = 0.0, 1.0
    for _ in range(14):
        mid = (lo + hi) / 2
        if winding_number(lerp(p0, c, mid), a.polys) != 0:
            lo = mid
        else:
            hi = mid
    side[end_i] = lerp(p0, c, lo)


def _build_ribbons(a: "Analysis", strokes: List["Stroke"], corners, factor: float, step: float,
                   abs_width: Optional[float] = None, fit_ink: bool = True,
                   include_dots: bool = True) -> List[List[Pt]]:
    """«Ленты» штрихов вокруг осей с толщиной factor * измеренная (или постоянной abs_width).
    Объединение лент даёт предсказанный контур. Для замкнутого штриха внешний контур идёт против
    часовой, внутренний по часовой, поэтому дыра остаётся дырой. fit_ink: торцы режутся по
    исходному контуру, торчащие углы подрезаются (для осей постоянной ширины это отключается)."""
    polys: List[List[Pt]] = []
    ends: Dict[Tuple[int, int], Tuple[Pt, Pt]] = {}

    def wf(w: float) -> float:
        return abs_width if abs_width else factor * w

    for s in (a.strokes if include_dots else []):  # точки: тот же контур, уменьшенный вокруг центра
        if s.kind == "dot" and s.patches:
            pts = s.patches[0]
            cx = sum(q[0] for q in pts) / len(pts)
            cy = sum(q[1] for q in pts) / len(pts)
            polys.append(_orient([(cx + (q[0] - cx) * factor, cy + (q[1] - cy) * factor) for q in pts], True))
    for s in strokes:
        pts = _axis_samples(s, step)
        closed = s.closed and len(pts) > 3
        if closed and dist(pts[0], pts[-1]) < 1.0:
            pts = pts[:-1]
        n = len(pts)
        left: List[Pt] = []
        right: List[Pt] = []
        for k, q in enumerate(pts):
            a_ = pts[(k - 1) % n] if closed else pts[max(k - 1, 0)]
            b_ = pts[(k + 1) % n] if closed else pts[min(k + 1, n - 1)]
            t = vunit(vsub(b_, a_)) or (1.0, 0.0)
            hw = 0.5 * wf(_width_at(s, q))
            left.append((q[0] - t[1] * hw, q[1] + t[0] * hw))
            right.append((q[0] + t[1] * hw, q[1] - t[0] * hw))
        if closed:
            outer, inner = (left, right) if abs(sum(vcross(left[i], left[(i + 1) % n]) for i in range(n))) >= \
                abs(sum(vcross(right[i], right[(i + 1) % n]) for i in range(n))) else (right, left)
            polys.append(_orient(outer, True))
            polys.append(_orient(inner, False))
        else:
            for idx, (L_i, R_i, k0, k1) in ((0, (0, 0, 0, 1)), (1, (-1, -1, -1, -2))):
                cap = s.caps.get(idx) if fit_ink else None
                if cap is None:
                    continue
                H, dc = cap
                t_out = vunit(vsub(pts[k0], pts[k1]))  # касательная наружу у этого конца
                if t_out is None:
                    continue
                for side in (left, right):
                    den = vcross(t_out, dc)
                    if abs(den) < 0.25:
                        continue  # торец почти параллелен оси: оставляем перпендикулярный срез
                    lam = vcross(vsub(H, side[L_i]), dc) / den
                    lim = 1.5 * wf(_width_at(s, pts[k0])) + 1.0
                    if abs(lam) <= lim:
                        side[L_i] = (side[L_i][0] + t_out[0] * lam, side[L_i][1] + t_out[1] * lam)
            if fit_ink and len(left) >= 2:
                for side in (left, right):
                    _clamp_to_ink(a, side, 0, 1)
                    _clamp_to_ink(a, side, -1, -2)
            polys.append(_orient(left + right[::-1], True))
            ends[(id(s), 0)] = (left[0], right[0])
            ends[(id(s), 1)] = (left[-1], right[-1])
    # угловые стыки: косой срез недостающего угла (митра)
    for X, e1, e2 in corners:
        if (id(e1.s), e1.idx) not in ends or (id(e2.s), e2.idx) not in ends:
            continue
        d1 = (-e1.t[0], -e1.t[1])
        d2 = (-e2.t[0], -e2.t[1])
        den = vcross(d1, d2)
        if abs(den) < 1e-6:
            continue
        best = None
        for c1 in ends[(id(e1.s), e1.idx)]:
            for c2 in ends[(id(e2.s), e2.idx)]:
                g = vsub(c2, c1)
                s_, u_ = vcross(g, d2) / den, vcross(g, d1) / den
                if s_ <= 1e-6 and u_ <= 1e-6 and (best is None or -(s_ + u_) < best[0]):
                    best = (-(s_ + u_), c1, c2, s_, u_)
        if best:
            _, c1, c2, s_, u_ = best
            lim = 2.0 * wf(max(e1.w, e2.w))
            M = (c1[0] + d1[0] * s_, c1[1] + d1[1] * s_)
            if abs(s_) <= lim and abs(u_) <= lim and winding_number(M, a.polys) != 0:
                polys.append(_orient([c1, M, c2, X], True))  # косой срез внешнего угла (при исходной толщине: если он внутри контура)
            else:
                polys.append(_orient([c1, c2, X], True))
    return polys


def _score(a: "Analysis", polys: List[List[Pt]], mean_w: float) -> Dict[str, float]:
    """Насколько предсказанный контур (ленты при исходной толщине) совпадает с существующим."""
    tol = max(3.0, 0.12 * mean_w)
    closed = [p + [p[0]] for p in polys if len(p) >= 3]
    boxes = [(min(q[0] for q in p) - tol, min(q[1] for q in p) - tol,
              max(q[0] for q in p) + tol, max(q[1] for q in p) + tol) for p in closed]

    def near_ribbon(pt: Pt) -> bool:
        for p, bb in zip(closed, boxes):
            if bb[0] <= pt[0] <= bb[2] and bb[1] <= pt[1] <= bb[3] and dist_to_polyline(pt, p) <= tol:
                return True
        return False

    spacing = max(4.0, mean_w / 8.0)
    tot = ok = 0
    for e in a.edges:
        n = max(1, int(e.length / spacing))
        for i in range(n):
            tot += 1
            if near_ribbon(e.point((i + 0.5) / n)):
                ok += 1
    sp2 = max(8.0, mean_w / 4.0)
    tot2 = over = 0
    for p in closed:
        for q in _resample(p[:-1], sp2):
            tot2 += 1
            if winding_number(q, a.polys) == 0 and min(dist_to_polyline(q, e.poly) for e in a.edges) > tol:
                over += 1
    return {"cover": ok / tot if tot else 1.0, "over": over / tot2 if tot2 else 0.0}


def connect_axes(a: "Analysis", upm: int, prm: Params) -> None:
    """Соединяет оси штрихов в общую ось и строит предсказанный контур."""
    a.joints, a.end_marks, a.pred_orig, a.pred_thin, a.scores = [], [], [], [], {}
    a.cn_strokes, a.cn_corners, a.mean_w = [], [], 0.0
    cn = _Connector(a, upm, prm)
    dots = [s for s in a.strokes if s.kind == "dot"]
    if not cn.strokes and not dots:
        return
    if cn.strokes and prm.join_axes:
        cn.run()
    allw = [s.w_mean for s in cn.strokes + dots]
    mean_w = sum(allw) / len(allw)
    a.cn_strokes, a.cn_corners, a.mean_w = cn.strokes, cn.corners, mean_w
    step = max(3.0, mean_w / 6.0)  # шаг по оси от толщины штриха, а не от размера em
    a.pred_orig = _build_ribbons(a, cn.strokes, cn.corners, 1.0, step)
    a.pred_thin = _build_ribbons(a, cn.strokes, cn.corners, prm.thin_ratio, step)
    a.scores = _score(a, a.pred_orig, mean_w)


# ============================================================================
#  Часть 3г. Сохранение в OTF: оси и тонкое начертание
# ============================================================================
#
#  Оси можно сохранить двумя способами: как центральные линии (команды move/line/curve без
#  расширения) или как узкие залитые ленты. CFF/Type 2 не имеет семантики открытого контура:
#  при заливке под-контур автоматически замыкается. Центральный режим предназначен для импорта
#  путей в редактор шрифтов, а не для корректного отображения в наборе.
#  Тонкое начертание строится из лент штрихов. Сохранение поддерживается для OTF с контурами
#  CFF; advance-ширины и остальные таблицы шрифта не меняются.

from fontTools.pens.recordingPen import RecordingPen  # noqa: E402
from fontTools.pens.t2CharStringPen import T2CharStringPen  # noqa: E402


@dataclass
class ExportOptions:
    mode: str = "thin"  # 'thin' | 'axes'
    thin_mode: str = "ratio"  # 'ratio': доля от исходной толщины | 'absolute': одинаковая толщина штрихов
    thin_ratio: float = 0.4
    thin_abs: float = 40.0  # ед. шрифта (для thin_mode == 'absolute')
    axis_width: float = 10.0  # ширина линии оси в ед. шрифта (для mode == 'axes' и axis_style == 'ribbon')
    axis_style: str = "ribbon"  # 'ribbon': прежние узкие ленты | 'centerline': исходные центральные пути
    fit_tol: float = 0.8  # допуск подгонки кривых, ед.
    corner_deg: float = 35.0  # излом ломаной не меньше этого угла считается углом
    use_pathops: bool = True  # объединять наложения (если установлен skia-pathops)
    name_suffix: str = ""  # добавляется к названию шрифта
    keep_unprocessed: bool = True  # глифы без найденных штрихов остаются как в исходном шрифте


@dataclass
class ExportReport:
    total: int = 0
    done: int = 0
    skipped: List[str] = field(default_factory=list)  # нет штрихов: оставлены как есть
    failed: List[str] = field(default_factory=list)  # ошибка обработки: оставлены как есть
    warnings: List[str] = field(default_factory=list)
    cancelled: bool = False


def _chord_params(pts: List[Pt]) -> List[float]:
    u = [0.0]
    for i in range(1, len(pts)):
        u.append(u[-1] + dist(pts[i - 1], pts[i]))
    tot = u[-1] or 1.0
    return [x / tot for x in u]


def _gen_bezier(pts: List[Pt], u: List[float], t1: Pt, t2: Pt) -> Cubic:
    p0, p3 = pts[0], pts[-1]
    c00 = c01 = c11 = x0 = x1 = 0.0
    for ui, p in zip(u, pts):
        b0, b1, b2, b3 = (1 - ui) ** 3, 3 * ui * (1 - ui) ** 2, 3 * ui * ui * (1 - ui), ui ** 3
        a1 = (t1[0] * b1, t1[1] * b1)
        a2 = (t2[0] * b2, t2[1] * b2)
        c00 += vdot(a1, a1)
        c01 += vdot(a1, a2)
        c11 += vdot(a2, a2)
        tmp = (p[0] - (p0[0] * (b0 + b1) + p3[0] * (b2 + b3)), p[1] - (p0[1] * (b0 + b1) + p3[1] * (b2 + b3)))
        x0 += vdot(a1, tmp)
        x1 += vdot(a2, tmp)
    det = c00 * c11 - c01 * c01
    al = (x0 * c11 - x1 * c01) / det if abs(det) > 1e-12 else 0.0
    ar = (c00 * x1 - c01 * x0) / det if abs(det) > 1e-12 else 0.0
    seg = dist(p0, p3)
    if al < 1e-6 * seg or ar < 1e-6 * seg:  # вырожденное решение: ручки по трети хорды
        al = ar = seg / 3.0
    return (p0, (p0[0] + t1[0] * al, p0[1] + t1[1] * al), (p3[0] + t2[0] * ar, p3[1] + t2[1] * ar), p3)


def _max_err(pts: List[Pt], bez: Cubic, u: List[float]) -> Tuple[float, int]:
    mx, split = 0.0, len(pts) // 2
    for i in range(1, len(pts) - 1):
        d = dist(bez_pt(bez, u[i]), pts[i])
        if d >= mx:
            mx, split = d, i
    return mx, split


def _reparam(pts: List[Pt], bez: Cubic, u: List[float]) -> List[float]:
    res = []
    for p, ui in zip(pts, u):
        q = bez_pt(bez, ui)
        d1 = bez_deriv(bez, ui)
        P0, P1, P2, P3 = bez
        d2 = (6 * ((1 - ui) * (P2[0] - 2 * P1[0] + P0[0]) + ui * (P3[0] - 2 * P2[0] + P1[0])),
              6 * ((1 - ui) * (P2[1] - 2 * P1[1] + P0[1]) + ui * (P3[1] - 2 * P2[1] + P1[1])))
        num = (q[0] - p[0]) * d1[0] + (q[1] - p[1]) * d1[1]
        den = d1[0] ** 2 + d1[1] ** 2 + (q[0] - p[0]) * d2[0] + (q[1] - p[1]) * d2[1]
        res.append(ui if abs(den) < 1e-12 else max(0.0, min(1.0, ui - num / den)))
    return res


def _fit_cubic(pts: List[Pt], t1: Pt, t2: Pt, err: float, out: List[Cubic], depth: int = 0) -> None:
    """Подгонка кубик Безье к точкам (алгоритм Шнайдера). t1, t2: единичные касательные у концов:
    t1 вперёд по кривой, t2 назад, к кривой."""
    n = len(pts)
    if n == 2:
        d = dist(pts[0], pts[1]) / 3.0
        out.append((pts[0], (pts[0][0] + t1[0] * d, pts[0][1] + t1[1] * d),
                    (pts[1][0] + t2[0] * d, pts[1][1] + t2[1] * d), pts[1]))
        return
    u = _chord_params(pts)
    bez = _gen_bezier(pts, u, t1, t2)
    emax, split = _max_err(pts, bez, u)
    if emax < err:
        out.append(bez)
        return
    if emax < err * 8:
        for _ in range(5):
            u = _reparam(pts, bez, u)
            bez = _gen_bezier(pts, u, t1, t2)
            emax, split = _max_err(pts, bez, u)
            if emax < err:
                out.append(bez)
                return
    if depth > 30 or n < 4:
        out.append(bez)
        return
    split = max(1, min(n - 2, split))
    ct = vunit(vsub(pts[split - 1], pts[split + 1])) or t1
    _fit_cubic(pts[:split + 1], t1, ct, err, out, depth + 1)
    _fit_cubic(pts[split:], (-ct[0], -ct[1]), t2, err, out, depth + 1)


def _run_to_segments(run: List[Pt], tol: float, t1: Optional[Pt] = None, t2: Optional[Pt] = None) -> List[Segment]:
    """Участок между двумя углами: прямая, если точки лежат на хорде, иначе кубики."""
    if len(run) == 2:
        return [("L", (run[0], run[1]))]
    chord = vsub(run[-1], run[0])
    L = math.hypot(chord[0], chord[1])
    if L > 1e-9 and all(abs(vcross(chord, vsub(q, run[0]))) / L <= tol for q in run[1:-1]) \
            and vdot(vsub(run[len(run) // 2], run[0]), chord) > 0:
        return [("L", (run[0], run[-1]))]
    t1 = t1 or vunit(vsub(run[1], run[0])) or (1.0, 0.0)
    t2 = t2 or vunit(vsub(run[-2], run[-1])) or (-1.0, 0.0)
    cubs: List[Cubic] = []
    _fit_cubic(run, t1, t2, tol, cubs)
    return [("C", c) for c in cubs]


def fit_polygon(pts: List[Pt], tol: float, corner_deg: float) -> List[Segment]:
    """Замкнутая ломаная в контур из отрезков и кубик: углы остаются углами."""
    p: List[Pt] = [pts[0]]
    for q in pts[1:]:
        if dist(q, p[-1]) > 1e-6:
            p.append(q)
    if len(p) > 1 and dist(p[0], p[-1]) < 1e-6:
        p.pop()
    n = len(p)
    if n < 3:
        return []

    def turn(i: int) -> float:
        a_, b_ = vsub(p[i], p[i - 1]), vsub(p[(i + 1) % n], p[i])
        return abs(angle_deg(a_, b_)) if vunit(a_) and vunit(b_) else 0.0

    corners = [i for i in range(n) if turn(i) >= corner_deg]
    segs: List[Segment] = []
    if not corners:  # гладкое кольцо: начало берём в точке 0, касательная по соседям
        t1 = vunit(vsub(p[1], p[-1])) or (1.0, 0.0)
        return _run_to_segments(p + [p[0]], tol, t1, (-t1[0], -t1[1]))
    k0 = corners[0]
    p = p[k0:] + p[:k0]
    cs = sorted((c - k0) % n for c in corners)
    bounds = cs + [n]
    for a_, b_ in zip(bounds, bounds[1:]):
        run = p[a_:b_ + 1] if b_ < n else p[a_:] + [p[0]]
        if len(run) >= 2:
            segs.extend(_run_to_segments(run, tol))
    return segs


def _sample_curve(op: str, args: tuple, start: Pt) -> List[Pt]:
    if op == "curveTo":
        c = (start, args[0], args[1], args[2])
        return flatten(c, 8)[1:]
    if op == "qCurveTo":  # кривые TrueType: цепочка внеконтурных точек, конец последний
        pts = [start] + list(args)
        out: List[Pt] = []
        cur = start
        offs = list(args[:-1])
        end = args[-1]
        for i, off in enumerate(offs):
            nxt = end if i == len(offs) - 1 else lerp(off, offs[i + 1], 0.5)
            for k in range(1, 7):
                t = k / 6
                out.append(((1 - t) ** 2 * cur[0] + 2 * (1 - t) * t * off[0] + t * t * nxt[0],
                            (1 - t) ** 2 * cur[1] + 2 * (1 - t) * t * off[1] + t * t * nxt[1]))
            cur = nxt
        return out
    return [args[-1]]


def union_polygons(polys: List[List[Pt]]) -> List[List[Pt]]:
    """Объединение лент по правилу ненулевой намотки (skia-pathops). Внешние контуры идут против
    часовой стрелки, дыры по часовой."""
    path = pathops.Path()
    pen = path.getPen()
    for poly in polys:
        pen.moveTo(poly[0])
        for q in poly[1:]:
            pen.lineTo(q)
        pen.closePath()
    path.simplify(fix_winding=True, keep_starting_points=False)
    rec = RecordingPen()
    path.draw(rec)
    out: List[List[Pt]] = []
    cur: List[Pt] = []
    for op, args in rec.value:
        if op == "moveTo":
            cur = [args[0]]
        elif op == "lineTo":
            cur.append(args[0])
        elif op in ("curveTo", "qCurveTo") and cur:
            cur.extend(_sample_curve(op, args, cur[-1]))
        elif op in ("closePath", "endPath"):
            if len(cur) > 2:
                out.append(cur)
            cur = []
    return out


def circle_segments(cx: float, cy: float, r: float) -> List[Segment]:
    k = 0.5522847498 * r
    return [
        ("C", ((cx + r, cy), (cx + r, cy + k), (cx + k, cy + r), (cx, cy + r))),
        ("C", ((cx, cy + r), (cx - k, cy + r), (cx - r, cy + k), (cx - r, cy))),
        ("C", ((cx - r, cy), (cx - r, cy - k), (cx - k, cy - r), (cx, cy - r))),
        ("C", ((cx, cy - r), (cx + k, cy - r), (cx + r, cy - k), (cx + r, cy))),
    ]


def glyph_centerline_paths(a: "Analysis") -> List[List[Segment]]:
    """Возвращает найденные оси как отдельные открытые пути из L/C-сегментов.

    Геометрия берётся непосредственно из a.cn_strokes[*].axis: без расширения,
    смещения, построения лент, объединения pathops и повторной подгонки полигонов.
    Разрывы между соседними кубиками не соединяются искусственно.
    """
    paths: List[List[Segment]] = []
    for stroke in a.cn_strokes:
        chain: List[Segment] = []
        last: Optional[Pt] = None
        for cubic in stroke.axis:
            if last is not None and dist(last, cubic[0]) > 0.01:
                if chain:
                    paths.append(chain)
                chain = []
            # У прямых осей контрольные точки лежат на хорде; записываем их как L,
            # чтобы не раздувать количество команд в CFF.
            if is_flat(cubic, 1e-4):
                seg: Segment = ("L", (cubic[0], cubic[3]))
            else:
                seg = ("C", cubic)
            if not chain:
                chain = [seg]
            else:
                chain.append(seg)
            last = cubic[3]
        if chain:
            paths.append(chain)
    return paths


def _centerline_path_is_straight(path: List[Segment], tol: float = 1e-4) -> bool:
    """Проверяет, что весь путь состоит из прямых на одной общей линии."""
    if not path or any(kind != "L" for kind, _ in path):
        return False
    points = [path[0][1][0]] + [pts[1] for _kind, pts in path]
    if len(points) < 3:
        return True
    chord = vsub(points[-1], points[0])
    length = math.hypot(chord[0], chord[1])
    if length < 1e-9:
        return all(dist(points[0], q) <= tol for q in points[1:])
    return all(abs(vcross(chord, vsub(q, points[0]))) / length <= tol for q in points[1:-1])


def centerlines_to_charstring(paths: List[List[Segment]], width: float, glyph_set, private, global_subrs):
    """Пишет оси в CFF CharString без построения лент.

    В Type 2 нет семантики открытого контура: каждый subpath автоматически замыкается.
    Прямые оси оставляются одиночными путями: их автоматическое замыкание совпадает с той же
    линией. Для кривых и ломаных путь возвращается по себе в обратном направлении, чтобы
    закрывающей хорды и залитой области не появлялось. Поэтому у криволинейных осей в редакторе
    будут совпадающие проходы вперёд и назад — ограничение CFF/OTF без открытых контуров.
    """
    pen = T2CharStringPen(width, glyph_set, roundTolerance=0.0)
    for path in paths:
        if not path:
            continue
        first = path[0][1][0]
        pen.moveTo(first)
        for kind, pts in path:
            if kind == "L":
                pen.lineTo(pts[1])
            else:
                pen.curveTo(pts[1], pts[2], pts[3])
        if not _centerline_path_is_straight(path):
            # Возврат точно по той же геометрии: не появляется закрывающая хорда.
            for kind, pts in reversed(path):
                if kind == "L":
                    pen.lineTo(pts[0])
                else:
                    pen.curveTo(pts[2], pts[1], pts[0])
    return pen.getCharString(private, global_subrs)


def glyph_contours(a: "Analysis", opts: ExportOptions, warn: Optional[List[str]] = None) -> List[List[Segment]]:
    """Контуры (отрезки и кубики) для сохранения: оси или тонкое начертание этого глифа."""
    dots = [s for s in a.strokes if s.kind == "dot"]
    if not a.cn_strokes and not dots:
        return []
    mean_w = a.mean_w or 100.0
    contours: List[List[Segment]] = []
    if opts.mode == "axes":
        polys = _build_ribbons(a, a.cn_strokes, [], 1.0, max(2.0, opts.axis_width), opts.axis_width,
                               fit_ink=False, include_dots=False)
        for d in dots:  # точка: маленький круг в её центре
            pts = d.patches[0]
            contours.append(circle_segments(sum(q[0] for q in pts) / len(pts), sum(q[1] for q in pts) / len(pts),
                                            opts.axis_width))
    else:
        absw = opts.thin_abs if opts.thin_mode == "absolute" else None
        factor = opts.thin_ratio if absw is None else 1.0
        polys = _build_ribbons(a, a.cn_strokes, a.cn_corners, factor, max(2.0, mean_w / 14.0), absw,
                               fit_ink=True, include_dots=False)
        for d in dots:  # точка: её контур, уменьшенный вокруг центра
            if d.contour < 0:
                continue
            pts = d.patches[0]
            cx, cy = sum(q[0] for q in pts) / len(pts), sum(q[1] for q in pts) / len(pts)
            k = factor if absw is None else max(0.05, min(1.0, 1.1 * absw / max(d.w_mean, 1.0)))
            sc = lambda q: (cx + (q[0] - cx) * k, cy + (q[1] - cy) * k)
            contours.append([("L", (sc(e.pts[0]), sc(e.pts[3]))) if e.kind == "line"
                             else ("C", tuple(sc(q) for q in e.pts)) for e in a.contours[d.contour]])
    polys = [pl for pl in polys if len(pl) >= 3]
    if opts.use_pathops and HAVE_PATHOPS and polys:
        try:
            polys = union_polygons(polys)
        except Exception as ex:  # без объединения: наложения остаются отдельными контурами
            if warn is not None:
                warn.append(f"объединение наложений не удалось: {ex}")
    for pl in polys:
        segs = fit_polygon(pl, opts.fit_tol, opts.corner_deg)
        if segs:
            contours.append(segs)
    return contours


def contours_to_charstring(contours: List[List[Segment]], width: float, glyph_set, private, global_subrs):
    pen = T2CharStringPen(width, glyph_set)
    for c in contours:
        pen.moveTo(c[0][1][0])
        for kind, pts in c:
            if kind == "L":
                pen.lineTo(pts[1])
            else:
                pen.curveTo(pts[1], pts[2], pts[3])
        pen.closePath()
    return pen.getCharString(private, global_subrs)


def _apply_name_suffix(font: TTFont, suffix: str) -> None:
    suffix = suffix.strip()
    if not suffix:
        return
    ps_suffix = "".join(ch for ch in suffix if ch.isalnum())
    for rec in list(font["name"].names):
        if rec.nameID in (1, 4, 16):
            font["name"].setName(rec.toUnicode() + " " + suffix, rec.nameID, rec.platformID, rec.platEncID, rec.langID)
        elif rec.nameID in (3, 6):
            font["name"].setName(rec.toUnicode() + ps_suffix, rec.nameID, rec.platformID, rec.platEncID, rec.langID)
    try:
        cff = font["CFF "].cff
        cff.fontNames[0] = cff.fontNames[0] + ps_suffix
        top = cff.topDictIndex[0]
        for attr in ("FullName", "FamilyName"):
            if hasattr(top, attr):
                setattr(top, attr, getattr(top, attr) + " " + suffix)
    except Exception:
        pass


def export_font(src_path: str, out_path: str, profiles, opts: ExportOptions, progress=None) -> ExportReport:
    """Сохраняет новый OTF с центральными линиями, узкими лентами или тонким начертанием.
    progress(i, n, имя_глифа) -> bool: False отменяет работу."""
    if os.path.abspath(src_path) == os.path.abspath(out_path):
        raise ValueError("Нельзя перезаписывать исходный шрифт. Выберите другое имя файла.")
    font = TTFont(src_path)  # не ленивая загрузка: заменяем контуры в таблице CFF
    if "CFF " not in font:
        raise ValueError("Сохранение поддержано для OTF с контурами CFF; у этого шрифта их нет (TrueType или CFF2).")
    store = FontStore(src_path)
    top = font["CFF "].cff.topDictIndex[0]
    cstrings = top.CharStrings
    gsubrs = font["CFF "].cff.GlobalSubrs
    order = font.getGlyphOrder()
    rep = ExportReport(total=len(order))
    prof = profiles if isinstance(profiles, dict) else {"straight": profiles, "curved": profiles}
    centerline_mode = opts.mode == "axes" and opts.axis_style == "centerline"
    for gid, name in enumerate(order):
        if progress is not None and not progress(gid, len(order), name):
            rep.cancelled = True
            break
        try:
            raw = store.raw(name, prof["curved"].remove_overlaps)
            if not raw.contours:
                continue
            a = analyze(raw, store.upm, prof)
            if centerline_mode:
                paths = glyph_centerline_paths(a)
                if not paths:
                    rep.skipped.append(name)
                    if opts.keep_unprocessed:
                        continue
            else:
                warns: List[str] = []
                contours = glyph_contours(a, opts, warns)
                rep.warnings.extend(f"{name}: {w}" for w in warns[:1])
                if not contours:
                    rep.skipped.append(name)
                    if opts.keep_unprocessed:
                        continue
            if hasattr(top, "FDSelect") and getattr(top, "FDSelect", None) is not None and hasattr(top, "FDArray"):
                private = top.FDArray[top.FDSelect[gid]].Private
            else:
                private = top.Private
            if centerline_mode:
                cstrings[name] = centerlines_to_charstring(paths, raw.advance, store.glyph_set, private, gsubrs)
            else:
                cstrings[name] = contours_to_charstring(contours, raw.advance, store.glyph_set, private, gsubrs)
            rep.done += 1
        except Exception as ex:
            rep.failed.append(f"{name}: {type(ex).__name__}: {ex}")
    if not rep.cancelled:
        _apply_name_suffix(font, opts.name_suffix)
        font.save(out_path)
    return rep


def _svg_path_data(path: List[Segment]) -> str:
    """SVG d для одного сегментированного пути в координатах шрифта."""
    if not path:
        return ""
    p0 = path[0][1][0]
    parts = [f"M {p0[0]:.4f} {p0[1]:.4f}"]
    for kind, pts in path:
        if kind == "L":
            p = pts[1]
            parts.append(f"L {p[0]:.4f} {p[1]:.4f}")
        else:
            p1, p2, p3 = pts[1], pts[2], pts[3]
            parts.append(f"C {p1[0]:.4f} {p1[1]:.4f}, {p2[0]:.4f} {p2[1]:.4f}, {p3[0]:.4f} {p3[1]:.4f}")
    return " ".join(parts)


def export_svg(src_path: str, out_path: str, profiles, axis_width: float = 12.0, progress=None) -> ExportReport:
    """Экспортирует все найденные оси в один SVG-лист: отдельная группа на каждый глиф.
    Пути остаются линиями SVG (fill=none, stroke задан только для просмотра).
    """
    from xml.sax.saxutils import escape, quoteattr

    if os.path.abspath(src_path) == os.path.abspath(out_path):
        raise ValueError("Нельзя перезаписывать исходный шрифт. Выберите другое имя файла.")
    store = FontStore(src_path)
    order = store.display_order
    prof = profiles if isinstance(profiles, dict) else {"straight": profiles, "curved": profiles}
    rep = ExportReport(total=len(order))
    pad = store.upm * 0.18
    max_advance = max([store.raw(n, prof["curved"].remove_overlaps).advance for n in order] or [store.upm])
    max_ink_width = store.upm
    for name in order:
        bb = store.raw(name, prof["curved"].remove_overlaps).bbox()
        if bb:
            max_ink_width = max(max_ink_width, bb[2] - bb[0])
    cell_w = max(max_advance, max_ink_width) + 2 * pad
    cell_h = max(store.ascent - store.descent, store.upm) + 2 * pad
    cols = max(1, int(math.ceil(math.sqrt(max(1, len(order))))))
    rows = max(1, int(math.ceil(len(order) / cols)))
    view_w, view_h = cols * cell_w, rows * cell_h
    sw = max(1.0, float(axis_width))
    groups: List[str] = []
    for i, name in enumerate(order):
        if progress is not None and not progress(i, len(order), name):
            rep.cancelled = True
            break
        try:
            raw = store.raw(name, prof["curved"].remove_overlaps)
            if not raw.contours:
                rep.skipped.append(name)
                continue
            a = analyze(raw, store.upm, prof)
            paths = glyph_centerline_paths(a)
            if not paths:
                rep.skipped.append(name)
                continue
            bb = raw.bbox() or (0.0, store.descent, raw.advance, store.ascent)
            col, row = i % cols, i // cols
            tx = col * cell_w + pad - bb[0]
            baseline = row * cell_h + pad + store.ascent
            cp = store.cp_of.get(name)
            cp_attr = f' data-unicode="U+{cp:04X}"' if cp is not None else ""
            # transform сохраняет математическую систему координат глифа (Y вверх).
            paths_svg = []
            for path in paths:
                d = _svg_path_data(path)
                if d:
                    paths_svg.append(f'<path d={quoteattr(d)} />')
            if paths_svg:
                groups.append(
                    f'<g id="glyph-{i:05d}" data-glyph-name={quoteattr(name)}{cp_attr}>'
                    f'<title>{escape(name)}</title>'
                    f'<g transform="translate({tx:.4f} {baseline:.4f}) scale(1 -1)">'
                    + "".join(paths_svg) + "</g></g>"
                )
            rep.done += 1
        except Exception as ex:
            rep.failed.append(f"{name}: {type(ex).__name__}: {ex}")
    if not rep.cancelled:
        header = ('<?xml version="1.0" encoding="UTF-8"?>\n'
                  f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {view_w:.4f} {view_h:.4f}" '
                  f'width="{view_w:.4f}" height="{view_h:.4f}">\n'
                  f'<g fill="none" stroke="#171717" stroke-width="{sw:.4f}" '
                  'stroke-linecap="round" stroke-linejoin="round">\n')
        with open(out_path, "w", encoding="utf-8", newline="\n") as f:
            f.write(header)
            f.write("\n".join(groups))
            f.write("\n</g>\n</svg>\n")
    return rep


# ============================================================================
#  Часть 3д. Экспорт настоящих открытых осей: UFO, Designspace, Glyphs
# ============================================================================


def _font_name(font: TTFont, name_id: int, fallback: str = "") -> str:
    try:
        value = font["name"].getDebugName(name_id)
        return value.strip() if value else fallback
    except Exception:
        return fallback


def _make_ufo_info(store: FontStore):
    """Минимальные метаданные UFO из исходного шрифта.

    Оригинальные контуры и признаки шрифта не копируются: UFO содержит только
    открытые центральные линии, но сохраняет имена, кодировки, метрики и UPM.
    """
    font = store.font
    family = _font_name(font, 16, _font_name(font, 1, store.title or "StrokeLab"))
    style = "Centerlines"
    os2 = font["OS/2"] if "OS/2" in font else None
    post = font["post"] if "post" in font else None
    head = font["head"] if "head" in font else None
    revision = float(getattr(head, "fontRevision", 1.0) or 1.0)
    version_major = max(1, int(revision))
    version_minor = max(0, min(999, int(round((revision - int(revision)) * 1000))))

    values = {
        "familyName": family,
        "styleName": style,
        "unitsPerEm": int(store.upm),
        "ascender": float(store.ascent),
        "descender": float(store.descent),
        "versionMajor": version_major,
        "versionMinor": version_minor,
        "postscriptFontName": _font_name(font, 6, ""),
        "postscriptFullName": f"{family} {style}",
        "openTypeNamePreferredFamilyName": family,
        "openTypeNamePreferredSubfamilyName": style,
        "openTypeNameVersion": _font_name(font, 5, f"Version {version_major}.{version_minor:03d}"),
        "copyright": _font_name(font, 0, ""),
        "trademark": _font_name(font, 7, ""),
        "openTypeNameDesigner": _font_name(font, 9, ""),
        "openTypeNameManufacturer": _font_name(font, 8, ""),
        "openTypeNameLicense": _font_name(font, 13, ""),
        "openTypeNameLicenseURL": _font_name(font, 14, ""),
        "italicAngle": float(getattr(post, "italicAngle", 0.0) or 0.0),
        "openTypeOS2WeightClass": int(getattr(os2, "usWeightClass", 400) or 400),
        "openTypeOS2WidthClass": int(getattr(os2, "usWidthClass", 5) or 5),
    }
    if os2 is not None:
        xh = int(getattr(os2, "sxHeight", 0) or 0)
        ch = int(getattr(os2, "sCapHeight", 0) or 0)
        if xh > 0:
            values["xHeight"] = xh
        if ch > 0:
            values["capHeight"] = ch
    # UFOWriter принимает объект с атрибутами и сам отбирает допустимые поля.
    return SimpleNamespace(**{k: v for k, v in values.items() if v not in (None, "")})


def _draw_centerline_paths_to_point_pen(paths: List[List[Segment]], point_pen) -> None:
    """Передаёт оси в UFO PointPen как открытые пути GLIF.

    В отличие от CFF/OTF, GLIF прямо поддерживает beginPath/endPath без
    closePath. Для каждой центральной линии создаётся ровно один путь.
    """
    for path in paths:
        if not path:
            continue
        point_pen.beginPath()
        start = path[0][1][0]
        point_pen.addPoint(start, segmentType="move")
        for kind, pts in path:
            if kind == "L":
                point_pen.addPoint(pts[1], segmentType="line")
            else:
                point_pen.addPoint(pts[1], segmentType=None)
                point_pen.addPoint(pts[2], segmentType=None)
                point_pen.addPoint(pts[3], segmentType="curve")
        point_pen.endPath()


def _unicode_map_for_store(store: FontStore) -> Dict[str, List[int]]:
    reverse: Dict[str, List[int]] = {}
    for cp, name in (store.font.getBestCmap() or {}).items():
        reverse.setdefault(name, []).append(int(cp))
    for values in reverse.values():
        values.sort()
    return reverse


def export_ufo(src_path: str, out_path: str, profiles, progress=None) -> ExportReport:
    """Экспортирует все глифы в UFO 3 с настоящими открытыми центральными линиями.

    Глифы без найденных осей сохраняются пустыми, но с исходной advance width и
    Unicode. Это сохраняет сетку глифов и кодировку, не подмешивая исходные заливки.
    """
    from fontTools.ufoLib import UFOWriter

    if os.path.abspath(src_path) == os.path.abspath(out_path):
        raise ValueError("Нельзя перезаписывать исходный шрифт.")
    if not out_path.lower().endswith(".ufo"):
        raise ValueError("Для UFO укажите имя с расширением .ufo.")
    if os.path.exists(out_path):
        raise FileExistsError(f"Папка или файл уже существует: {out_path}\nВыберите новое имя, чтобы не перезаписать данные.")
    parent = os.path.dirname(os.path.abspath(out_path))
    if not os.path.isdir(parent):
        raise FileNotFoundError(f"Папка назначения не существует: {parent}")

    store = FontStore(src_path)
    order = store.font.getGlyphOrder()
    prof = profiles if isinstance(profiles, dict) else {"straight": profiles, "curved": profiles}
    unicode_map = _unicode_map_for_store(store)
    rep = ExportReport(total=len(order))

    temp_root = tempfile.mkdtemp(prefix=".strokelab_ufo_", dir=parent)
    staged_path = os.path.join(temp_root, os.path.basename(out_path))
    writer = None
    try:
        writer = UFOWriter(staged_path, formatVersion=3)
        glyph_set = writer.getGlyphSet()
        writer.writeInfo(_make_ufo_info(store))
        writer.writeLib({
            "public.glyphOrder": list(order),
            "com.strokelab.centerlineExport": {
                "sourceFont": os.path.basename(src_path),
                "description": "Open centerline paths exported by StrokeLab",
                "pathSemantics": "open",
            },
        })

        for i, name in enumerate(order):
            if progress is not None and not progress(i, len(order), name):
                rep.cancelled = True
                break
            glyph_obj = store.glyph_set[name]
            width = float(getattr(glyph_obj, "width", 0) or 0)
            paths: List[List[Segment]] = []
            try:
                raw = store.raw(name, prof["curved"].remove_overlaps)
                if raw.contours:
                    analysis = analyze(raw, store.upm, prof)
                    paths = glyph_centerline_paths(analysis)
            except Exception as ex:
                rep.failed.append(f"{name}: {type(ex).__name__}: {ex}")

            glyph_data = SimpleNamespace(
                width=width,
                unicodes=unicode_map.get(name, []),
                note="StrokeLab centerlines; open paths" if paths else "No centerline detected",
            )
            draw_func = (lambda point_pen, p=paths: _draw_centerline_paths_to_point_pen(p, point_pen)) if paths else None
            try:
                glyph_set.writeGlyph(name, glyph_data, drawPointsFunc=draw_func)
                if paths:
                    rep.done += 1
                else:
                    rep.skipped.append(name)
            except Exception as ex:
                rep.failed.append(f"{name}: запись UFO: {type(ex).__name__}: {ex}")

        if not rep.cancelled:
            glyph_set.writeContents()
            writer.writeLayerContents()
        writer.close()
        writer = None

        if not rep.cancelled:
            os.rename(staged_path, out_path)
    except Exception:
        if writer is not None:
            try:
                writer.close()
            except Exception:
                pass
        raise
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)
    return rep


def export_designspace(src_path: str, out_path: str, profiles, progress=None) -> ExportReport:
    """Экспортирует DesignSpace 5.0 и связанный с ним один UFO-мастер.

    Это корректная одностилевая оболочка Designspace+UFO. Она не создаёт
    интерполяцию: для переменного шрифта потребуются отдельные UFO для мастеров.
    """
    from fontTools.designspaceLib import DesignSpaceDocument, SourceDescriptor

    if not out_path.lower().endswith(".designspace"):
        raise ValueError("Для Designspace укажите расширение .designspace.")
    if os.path.abspath(src_path) == os.path.abspath(out_path):
        raise ValueError("Путь назначения совпадает с исходным файлом.")
    final_ds = os.path.abspath(out_path)
    parent = os.path.dirname(final_ds)
    stem = os.path.splitext(os.path.basename(final_ds))[0]
    final_ufo = os.path.join(parent, stem + ".ufo")
    if os.path.exists(final_ds) or os.path.exists(final_ufo):
        raise FileExistsError(
            "Файл Designspace или его UFO-мастер уже существует. "
            "Выберите другое имя, чтобы не перезаписать данные."
        )
    if not os.path.isdir(parent):
        raise FileNotFoundError(f"Папка назначения не существует: {parent}")

    store = FontStore(src_path)
    family = _font_name(store.font, 16, _font_name(store.font, 1, store.title or "StrokeLab"))
    temp_root = tempfile.mkdtemp(prefix=".strokelab_designspace_", dir=parent)
    staged_ufo = os.path.join(temp_root, os.path.basename(final_ufo))
    staged_ds = os.path.join(temp_root, os.path.basename(final_ds))
    try:
        rep = export_ufo(src_path, staged_ufo, profiles, progress)
        if rep.cancelled:
            return rep

        doc = DesignSpaceDocument()
        doc.formatVersion = "5.0"
        source = SourceDescriptor()
        source.filename = os.path.basename(final_ufo)  # путь относительный к .designspace
        source.name = "Centerlines"
        source.familyName = family
        source.styleName = "Centerlines"
        source.location = {}
        source.copyInfo = True
        doc.addSource(source)
        doc.default = source.name
        doc.lib["com.strokelab.centerlineExport"] = {
            "description": "Single-source DesignSpace containing open centerline paths",
            "note": "Add additional master UFOs and axes before using for interpolation.",
        }
        doc.write(staged_ds)
        # Переносим оба результата как одну операцию: при сбое записи .designspace
        # не оставляем рядом UFO-мастер без соответствующего документа.
        ufo_moved = False
        try:
            os.rename(staged_ufo, final_ufo)
            ufo_moved = True
            os.rename(staged_ds, final_ds)
        except Exception:
            if ufo_moved and os.path.isdir(final_ufo):
                shutil.rmtree(final_ufo, ignore_errors=True)
            raise
        rep.warnings.append(
            "Создан один источник UFO без осей вариации. Для интерполяции добавьте другие мастера и оси."
        )
        return rep
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def export_glyphs(src_path: str, out_path: str, profiles, progress=None) -> ExportReport:
    """Экспортирует открытые оси в Glyphs 3 через промежуточный UFO 3.

    glyphsLib и ufoLib2 — необязательные зависимости. Переход UFO -> Glyphs
    сохраняет открытый характер путей и кубические сегменты.
    """
    if not out_path.lower().endswith(".glyphs"):
        raise ValueError("Для Glyphs укажите расширение .glyphs.")
    if os.path.abspath(src_path) == os.path.abspath(out_path):
        raise ValueError("Путь назначения совпадает с исходным файлом.")
    if os.path.exists(out_path):
        raise FileExistsError(f"Файл уже существует: {out_path}\nВыберите новое имя, чтобы не перезаписать данные.")
    try:
        import glyphsLib
        import ufoLib2
    except Exception as ex:
        raise RuntimeError(
            "Для экспорта Glyphs нужны glyphsLib и ufoLib2. Установите их командой:\n"
            "python -m pip install --upgrade glyphsLib ufoLib2"
        ) from ex

    parent = os.path.dirname(os.path.abspath(out_path))
    if not os.path.isdir(parent):
        raise FileNotFoundError(f"Папка назначения не существует: {parent}")
    temp_root = tempfile.mkdtemp(prefix=".strokelab_glyphs_", dir=parent)
    temp_ufo = os.path.join(temp_root, "Centerlines.ufo")
    temp_glyphs = os.path.join(temp_root, "Centerlines.glyphs")
    try:
        rep = export_ufo(src_path, temp_ufo, profiles, progress)
        if rep.cancelled:
            return rep
        ufo = ufoLib2.Font.open(temp_ufo)
        gsfont = glyphsLib.to_glyphs([ufo], ufo_module=ufoLib2)
        # Записываем современный текстовый формат Glyphs 3.
        gsfont.format_version = 3
        try:
            gsfont.note = "StrokeLab centerline source. Paths are intentionally open."
        except Exception:
            pass
        gsfont.save(temp_glyphs)
        os.rename(temp_glyphs, out_path)
        return rep
    except Exception as ex:
        if isinstance(ex, RuntimeError) and "glyphsLib" in str(ex):
            raise
        raise RuntimeError(f"Не удалось сформировать Glyphs-файл: {type(ex).__name__}: {ex}") from ex
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


# ============================================================================
#  Часть 4. Интерфейс (PySide6)
# ============================================================================

from PySide6.QtCore import QAbstractListModel, QModelIndex, QPointF, QRectF, QSize, Qt, QTimer, Signal  # noqa: E402
from PySide6.QtGui import (  # noqa: E402
    QColor,
    QFont,
    QIcon,
    QKeySequence,
    QPainter,
    QPainterPath,
    QPen,
    QPixmap,
    QPolygonF,
    QShortcut,
)
from PySide6.QtWidgets import (  # noqa: E402
    QApplication,
    QCheckBox,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QGridLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListView,
    QMainWindow,
    QMessageBox,
    QPlainTextEdit,
    QProgressDialog,
    QPushButton,
    QScrollArea,
    QSizePolicy,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

VERTEX_NAMES = {
    "convex": "выпуклый угол",
    "concave": "вогнутый угол",
    "smooth": "гладкий узел",
    "inflection": "перегиб",
}
STROKE_NAMES = {"line": "прямой штрих", "arc": "штрих-дуга",
                "blob": "пятно (не штрих: засечка, стык)", "dot": "точка (i, j, двоеточие)"}


def build_qpath(contours: List[List[Segment]], s: float, ox: float, oy: float) -> QPainterPath:
    """Контуры глифа в QPainterPath; координаты шрифта -> экран: (ox + x*s, oy - y*s)."""
    path = QPainterPath()
    path.setFillRule(Qt.FillRule.WindingFill)
    for c in contours:
        if not c:
            continue
        x0, y0 = c[0][1][0]
        path.moveTo(ox + x0 * s, oy - y0 * s)
        for kind, pts in c:
            if kind == "L":
                path.lineTo(ox + pts[1][0] * s, oy - pts[1][1] * s)
            else:
                path.cubicTo(
                    ox + pts[1][0] * s, oy - pts[1][1] * s,
                    ox + pts[2][0] * s, oy - pts[2][1] * s,
                    ox + pts[3][0] * s, oy - pts[3][1] * s,
                )
        path.closeSubpath()
    return path


class GlyphModel(QAbstractListModel):
    """Список глифов с миниатюрами, которые рисуются по требованию."""

    ICON = 56

    def __init__(self):
        super().__init__()
        self.store: Optional[FontStore] = None
        self.names: List[str] = []
        self._icons: Dict[str, QIcon] = {}

    def set_store(self, store: FontStore):
        self.beginResetModel()
        self.store = store
        self._icons.clear()
        self.names = list(store.display_order)
        self.endResetModel()

    def set_filter(self, text: str, by_name: bool = False):
        """Пустой фильтр: все глифы. По умолчанию ищет символы из строки
        (в порядке ввода). С by_name=True ищет глифы, в имени которых есть строка."""
        if not self.store:
            return
        self.beginResetModel()
        if not text:
            self.names = list(self.store.display_order)
        elif by_name:
            low = text.lower()
            self.names = [n for n in self.store.display_order if low in n.lower()]
        else:
            names: List[str] = []
            for ch in text:
                n = self.store.name_of_char.get(ch)
                if n and n not in names:
                    names.append(n)
            self.names = names
        self.endResetModel()

    def rowCount(self, parent=QModelIndex()):
        return 0 if parent.isValid() else len(self.names)

    def data(self, index, role=Qt.ItemDataRole.DisplayRole):
        if not index.isValid() or not self.store:
            return None
        name = self.names[index.row()]
        if role == Qt.ItemDataRole.DisplayRole:
            return name if len(name) <= 9 else name[:8] + "…"
        if role == Qt.ItemDataRole.DecorationRole:
            ic = self._icons.get(name)
            if ic is None:
                ic = self._make_icon(name)
                self._icons[name] = ic
            return ic
        if role == Qt.ItemDataRole.ToolTipRole:
            return self.store.describe(name)
        return None

    def _make_icon(self, name: str) -> QIcon:
        size = self.ICON
        try:
            raw = self.store.raw(name)
            pm = QPixmap(size, size)
            pm.fill(Qt.GlobalColor.transparent)
            p = QPainter(pm)
            p.setRenderHint(QPainter.RenderHint.Antialiasing, True)
            s = size / (self.store.upm * 1.25)
            ox = (size - raw.advance * s) / 2.0
            oy = size * 0.74
            p.fillPath(build_qpath(raw.contours, s, ox, oy), QColor(30, 30, 35))
            p.end()
            return QIcon(pm)
        except Exception:
            return QIcon()


class GlyphCanvas(QWidget):
    """Холст: заливка, пятна штрихов, оси, рёбра, нормали, узлы поверх глифа."""

    picked = Signal(object, object)  # (выбранное ребро или None, выбранный штрих или None)

    def __init__(self, interactive: bool = True):
        super().__init__()
        self.interactive = interactive  # миниатюры панели сравнения не реагируют на мышь
        self.setMinimumSize(520, 420) if interactive else self.setMinimumSize(0, 0)
        self.analysis: Optional[Analysis] = None
        self.ascent = 800.0
        self.descent = -200.0
        self.scale = 0.3
        self.origin = QPointF(100.0, 400.0)
        self.layers = {
            "fill": True, "strokes": True, "axes": True, "joints": True, "pred_thin": False, "clip_thin": True,
            "pred_orig": False, "edges": True, "leftover": True,
            "normals": False, "vertices": True, "labels": False, "metrics": True,
        }
        self.color_mode = "type"  # 'type' | 'index'
        self.sel_edge: Optional[Edge] = None
        self.sel_stroke: Optional[Stroke] = None
        self._press_pos: Optional[QPointF] = None
        self._last_pos: Optional[QPointF] = None
        self._dragging = False

    # --- данные --------------------------------------------------------------
    def set_analysis(self, a: Optional[Analysis], ascent: float, descent: float, refit: bool):
        self.analysis = a
        self.ascent, self.descent = ascent, descent
        self.sel_edge = None
        self.sel_stroke = None
        if refit:
            self.fit_view()
        self.update()

    def fit_view(self):
        a = self.analysis
        if a is None:
            return
        x0, y0, x1, y1 = a.bbox
        x0, x1 = min(x0, 0.0), max(x1, a.advance)
        y0, y1 = min(y0, 0.0), max(y1, 0.0)
        if y1 - y0 < 1 or x1 - x0 < 1:
            y0, y1 = self.descent, self.ascent
            x1 = max(x1, 100.0)
        margin = 48.0
        w, h = max(1.0, x1 - x0), max(1.0, y1 - y0)
        s = min((self.width() - 2 * margin) / w, (self.height() - 2 * margin) / h)
        self.scale = max(s, 1e-4)
        self.origin = QPointF(
            self.width() / 2.0 - (x0 + x1) / 2.0 * self.scale,
            self.height() / 2.0 + (y0 + y1) / 2.0 * self.scale,
        )
        self.update()

    # --- вспомогательное -----------------------------------------------------
    def tos(self, p: Pt) -> QPointF:
        return QPointF(self.origin.x() + p[0] * self.scale, self.origin.y() - p[1] * self.scale)

    def edge_color(self, e: Edge) -> QColor:
        if self.color_mode == "index":
            c = QColor.fromHsvF(((e.gid * 137.508) % 360.0) / 360.0, 0.78, 0.86)
        elif e.kind == "line":
            c = QColor(40, 120, 220)
        elif e.bulge == "convex":
            c = QColor(235, 120, 20)
        else:
            c = QColor(30, 160, 100)
        c.setAlpha(225)
        return c

    @staticmethod
    def stroke_color(s: Stroke, alpha: int) -> QColor:
        c = QColor(140, 140, 150) if s.kind == "blob" else QColor.fromHsvF(s.hue, 0.65, 0.95)
        c.setAlpha(alpha)
        return c

    @staticmethod
    def axis_color(s: Stroke) -> QColor:
        return QColor.fromHsvF(s.hue, 0.9, 0.42)

    def _edge_path(self, e: Edge) -> QPainterPath:
        path = QPainterPath(self.tos(e.pts[0]))
        if e.kind == "line":
            path.lineTo(self.tos(e.pts[3]))
        else:
            path.cubicTo(self.tos(e.pts[1]), self.tos(e.pts[2]), self.tos(e.pts[3]))
        return path

    def _poly_path(self, pts: List[Pt]) -> QPainterPath:
        path = QPainterPath(self.tos(pts[0]))
        for q in pts[1:]:
            path.lineTo(self.tos(q))
        path.closeSubpath()
        path.setFillRule(Qt.FillRule.WindingFill)
        return path

    def _axis_path(self, s: Stroke) -> QPainterPath:
        path = QPainterPath()
        for c in s.axis:
            path.moveTo(self.tos(c[0]))
            path.cubicTo(self.tos(c[1]), self.tos(c[2]), self.tos(c[3]))
        return path

    @staticmethod
    def _arrow(p: QPainter, a0: QPointF, a1: QPointF):
        p.drawLine(a0, a1)
        dx, dy = a1.x() - a0.x(), a1.y() - a0.y()
        l = math.hypot(dx, dy)
        if l < 1e-6:
            return
        ux, uy = -dx / l, -dy / l  # назад от острия
        for sgn in (1, -1):
            ang = math.radians(25 * sgn)
            c, s = math.cos(ang), math.sin(ang)
            hx, hy = ux * c - uy * s, ux * s + uy * c
            p.drawLine(a1, QPointF(a1.x() + hx * 6.5, a1.y() + hy * 6.5))

    # --- рисование -----------------------------------------------------------
    def paintEvent(self, _ev):
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        p.fillRect(self.rect(), QColor(252, 252, 250))
        a = self.analysis
        if a is None:
            p.setPen(QColor(130, 130, 130))
            p.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter,
                       "Откройте шрифт (Ctrl+O) или перетащите .otf в окно")
            return
        L = self.layers
        p.setBrush(Qt.BrushStyle.NoBrush)

        if L["metrics"]:
            pen = QPen(QColor(170, 170, 195))
            pen.setWidthF(1.0)
            pen.setStyle(Qt.PenStyle.DashLine)
            p.setPen(pen)
            oy = self.origin.y()
            for y in (0.0, self.ascent, self.descent):
                yy = oy - y * self.scale
                p.drawLine(QPointF(0, yy), QPointF(self.width(), yy))
            for x in (0.0, a.advance):
                xx = self.origin.x() + x * self.scale
                p.drawLine(QPointF(xx, 0), QPointF(xx, self.height()))

        if L["fill"]:
            path = build_qpath(a.raw_contours, self.scale, self.origin.x(), self.origin.y())
            p.fillPath(path, QColor(60, 60, 70, 40))
            pen = QPen(QColor(60, 60, 70, 120))
            pen.setWidthF(1.0)
            p.setPen(pen)
            p.drawPath(path)

        if L["pred_thin"] and a.pred_thin:  # предсказание тонкого начертания: объединение лент штрихов
            path = QPainterPath()
            path.setFillRule(Qt.FillRule.WindingFill)
            for poly in a.pred_thin:
                if len(poly) >= 3:
                    path.addPath(self._poly_path(poly))
            if L["clip_thin"]:  # тонкое начертание целиком лежит внутри исходного
                path = path.intersected(build_qpath(a.raw_contours, self.scale, self.origin.x(), self.origin.y()))
            p.fillPath(path, QColor(20, 90, 190, 125))

        if L["pred_orig"]:  # проверка: ленты при исходной толщине должны ложиться на контур
            pen = QPen(QColor(210, 0, 120, 210))
            pen.setWidthF(1.3)
            pen.setStyle(Qt.PenStyle.DashLine)
            p.setPen(pen)
            for poly in a.pred_orig:
                if len(poly) >= 3:
                    p.drawPath(self._poly_path(poly))

        if L["strokes"]:  # полупрозрачные пятна штрихов
            for s in a.strokes:
                sel = s is self.sel_stroke
                for patch in s.patches:
                    if len(patch) < 3:
                        continue
                    path = self._poly_path(patch)
                    p.fillPath(path, self.stroke_color(s, 150 if sel else 90))
                    pen = QPen(self.stroke_color(s, 230 if sel else 140))
                    pen.setWidthF(2.0 if sel else 1.0)
                    p.setPen(pen)
                    p.drawPath(path)

        if self.sel_edge is not None:
            pen = QPen(QColor(255, 205, 0, 170))
            pen.setWidthF(10.0)
            pen.setCapStyle(Qt.PenCapStyle.RoundCap)
            p.setPen(pen)
            p.drawPath(self._edge_path(self.sel_edge))

        if L["edges"]:
            for e in a.edges:
                pen = QPen(self.edge_color(e))
                pen.setWidthF(2.4 if L["strokes"] else 3.2)
                pen.setCapStyle(Qt.PenCapStyle.RoundCap)
                p.setPen(pen)
                p.drawPath(self._edge_path(e))

        if L["leftover"]:  # части рёбер, не вошедшие ни в один штрих
            pen = QPen(QColor(205, 0, 120, 210))
            pen.setWidthF(3.4)
            pen.setStyle(Qt.PenStyle.DotLine)
            p.setPen(pen)
            for e in a.edges:
                for t0, t1 in e.uncovered():
                    c = sub_cubic(e.pts, t0, t1)
                    path = QPainterPath(self.tos(c[0]))
                    path.cubicTo(self.tos(c[1]), self.tos(c[2]), self.tos(c[3]))
                    p.drawPath(path)

        if L["normals"]:
            pen = QPen(QColor(70, 70, 85, 210))
            pen.setWidthF(1.3)
            p.setPen(pen)
            for e in a.edges:
                ts = (0.5,) if (e.kind == "line" or e.length * self.scale < 40) else (0.2, 0.5, 0.8)
                for t in ts:
                    q = self.tos(e.point(t))
                    nx, ny = e.normal(t)
                    self._arrow(p, q, QPointF(q.x() + nx * 22, q.y() - ny * 22))

        if L["axes"]:  # оси штрихов; продления, мосты и достройки до торца пунктиром
            for s in a.strokes:
                if not s.axis:
                    continue
                col = self.axis_color(s)
                tags = s.tags or ["fit"] * len(s.axis)
                for c, tag in zip(s.axis, tags):
                    solid = tag in ("fit", "mid")  # mid: ось по середине контура в стыке
                    pen = QPen(col if solid else QColor(224, 128, 0))
                    pen.setWidthF(3.4 if s is self.sel_stroke else 2.6)
                    pen.setCapStyle(Qt.PenCapStyle.RoundCap)
                    if not solid:
                        pen.setStyle(Qt.PenStyle.DashLine)
                    p.setPen(pen)
                    path = QPainterPath(self.tos(c[0]))
                    path.cubicTo(self.tos(c[1]), self.tos(c[2]), self.tos(c[3]))
                    p.drawPath(path)
                if not s.closed:
                    p.setBrush(col)
                    p.setPen(Qt.PenStyle.NoPen)
                    for q in (s.axis[0][0], s.axis[-1][3]):
                        p.drawEllipse(self.tos(q), 3.6, 3.6)
                    p.setBrush(Qt.BrushStyle.NoBrush)

        if L["joints"]:  # соединения осей: T, L, Y, X, G1
            pen = QPen(QColor(30, 30, 35))
            pen.setWidthF(1.8)
            p.setPen(pen)
            p.setBrush(QColor(255, 255, 255, 235))
            for j in a.joints:
                c = self.tos(j.pos)
                r = 6.5
                if j.kind == "T":
                    p.drawPolygon(QPolygonF([QPointF(c.x() - r, c.y() - r * 0.8), QPointF(c.x() + r, c.y() - r * 0.8),
                                             QPointF(c.x(), c.y() + r)]))
                elif j.kind == "Y":
                    p.drawPolygon(QPolygonF([QPointF(c.x() - r, c.y() + r * 0.8), QPointF(c.x() + r, c.y() + r * 0.8),
                                             QPointF(c.x(), c.y() - r)]))
                elif j.kind == "L":
                    p.drawRect(QRectF(c.x() - r * 0.8, c.y() - r * 0.8, r * 1.6, r * 1.6))
                elif j.kind == "X":
                    p.drawEllipse(c, r, r)
                    p.drawLine(QPointF(c.x() - r * 0.7, c.y() - r * 0.7), QPointF(c.x() + r * 0.7, c.y() + r * 0.7))
                    p.drawLine(QPointF(c.x() - r * 0.7, c.y() + r * 0.7), QPointF(c.x() + r * 0.7, c.y() - r * 0.7))
                else:  # G1
                    p.drawEllipse(c, r * 0.8, r * 0.8)
            pen = QPen(QColor(235, 120, 0))
            pen.setWidthF(2.4)
            p.setPen(pen)
            p.setBrush(Qt.BrushStyle.NoBrush)
            for q, mk in a.end_marks:  # свободный конец: не торец и не соединён
                if mk == "open":
                    p.drawEllipse(self.tos(q), 8.5, 8.5)

        if L["vertices"]:
            white = QPen(QColor(255, 255, 255))
            white.setWidthF(1.5)
            p.setPen(white)
            for v in a.vertices:  # сначала мелкие гладкие узлы
                if v.kind != "smooth":
                    continue
                p.setBrush(QColor(120, 120, 130))
                p.drawEllipse(self.tos(v.pos), 2.6, 2.6)
            for v in a.vertices:
                c = self.tos(v.pos)
                if v.kind == "convex":
                    p.setBrush(QColor(225, 40, 40))
                    p.drawEllipse(c, 5.5, 5.5)
                elif v.kind == "concave":
                    p.setBrush(QColor(20, 20, 25))
                    p.drawRect(QRectF(c.x() - 4.5, c.y() - 4.5, 9, 9))
                elif v.kind == "inflection":
                    p.setBrush(QColor(150, 60, 200))
                    r = 6.0
                    p.drawPolygon(QPolygonF([
                        QPointF(c.x(), c.y() - r), QPointF(c.x() + r, c.y()),
                        QPointF(c.x(), c.y() + r), QPointF(c.x() - r, c.y()),
                    ]))
            p.setBrush(Qt.BrushStyle.NoBrush)

        if L["labels"]:
            f = QFont()
            f.setPointSizeF(8.5)
            f.setBold(True)
            p.setFont(f)
            for e in a.edges:  # номера рёбер
                c = self.tos(e.point(0.5))
                nx, ny = e.normal(0.5)
                x, y = c.x() - nx * 15, c.y() + ny * 15  # наружу от материала
                p.setPen(self.edge_color(e).darker(135))
                p.drawText(QRectF(x - 16, y - 9, 32, 18), Qt.AlignmentFlag.AlignCenter, str(e.gid))
            f.setPointSizeF(10.0)
            p.setFont(f)
            for j in a.joints:  # вид узла и опора в контуре (число вогнутых углов рядом)
                c = self.tos(j.pos)
                p.setPen(QColor(30, 30, 35))
                p.drawText(QRectF(c.x() + 7, c.y() - 22, 46, 16), Qt.AlignmentFlag.AlignLeft, f"{j.kind}{j.support}")
            for s in a.strokes:  # номера штрихов
                if s.axis:
                    c = self.tos(bez_pt(s.axis[len(s.axis) // 2], 0.5))
                    p.setPen(self.axis_color(s))
                else:
                    pts = s.patches[0]
                    c = self.tos((sum(q[0] for q in pts) / len(pts), sum(q[1] for q in pts) / len(pts)))
                    p.setPen(QColor(90, 90, 100))
                p.drawText(QRectF(c.x() - 18, c.y() - 22, 36, 18), Qt.AlignmentFlag.AlignCenter, f"S{s.sid}")

    # --- мышь ----------------------------------------------------------------
    def wheelEvent(self, ev):
        if not self.interactive:
            return
        dy = ev.angleDelta().y()
        if not dy:
            return
        f = 1.15 ** (dy / 120.0)
        f = max(0.002, min(300.0, self.scale * f)) / self.scale  # масштаб в пределах 0.002..300
        pos = ev.position()
        self.origin = QPointF(pos.x() + (self.origin.x() - pos.x()) * f,
                              pos.y() + (self.origin.y() - pos.y()) * f)
        self.scale *= f
        self.update()

    def mousePressEvent(self, ev):
        if not self.interactive:
            return
        self._press_pos = ev.position()
        self._last_pos = ev.position()
        self._dragging = False

    def mouseMoveEvent(self, ev):
        if not self.interactive or self._press_pos is None or self._last_pos is None:
            return
        pos = ev.position()
        if not self._dragging and math.hypot(pos.x() - self._press_pos.x(),
                                             pos.y() - self._press_pos.y()) > 4:
            self._dragging = True
        if self._dragging:
            self.origin = QPointF(self.origin.x() + pos.x() - self._last_pos.x(),
                                  self.origin.y() + pos.y() - self._last_pos.y())
            self.update()
        self._last_pos = pos

    def mouseReleaseEvent(self, ev):
        if not self.interactive:
            return
        if not self._dragging and ev.button() == Qt.MouseButton.LeftButton:
            self.sel_edge = self.pick_edge(ev.position())
            self.sel_stroke = self.pick_stroke(ev.position())
            self.picked.emit(self.sel_edge, self.sel_stroke)
            self.update()
        self._press_pos = None
        self._dragging = False

    def mouseDoubleClickEvent(self, _ev):
        if self.interactive:
            self.fit_view()

    def _to_font(self, pos: QPointF) -> Pt:
        return ((pos.x() - self.origin.x()) / self.scale, (self.origin.y() - pos.y()) / self.scale)

    def pick_edge(self, pos: QPointF) -> Optional[Edge]:
        if self.analysis is None or not self.layers["edges"]:
            return None
        pt = self._to_font(pos)
        best, best_d = None, 1e18
        for e in self.analysis.edges:
            d = dist_to_polyline(pt, e.poly)
            if d < best_d:
                best, best_d = e, d
        return best if best is not None and best_d * self.scale <= 10.0 else None

    def pick_stroke(self, pos: QPointF) -> Optional[Stroke]:
        a = self.analysis
        if a is None or not a.strokes or not (self.layers["strokes"] or self.layers["axes"]):
            return None
        pt = self._to_font(pos)
        best, best_area = None, 1e18
        for s in a.strokes:  # сначала: клик внутри пятна, выбираем самое маленькое
            for patch in s.patches:
                if len(patch) < 3 or winding_number(pt, [patch]) == 0:
                    continue
                xs, ys = [q[0] for q in patch], [q[1] for q in patch]
                area = (max(xs) - min(xs)) * (max(ys) - min(ys))
                if area < best_area:
                    best, best_area = s, area
        if best is not None:
            return best
        near, nd = None, 1e18  # иначе ближайшая ось
        for s in a.strokes:
            for c in s.axis:
                d = dist_to_polyline(pt, flatten(c, 12))
                if d < nd:
                    near, nd = s, d
        return near if near is not None and nd * self.scale <= 10.0 else None


class PreviewPanel(QWidget):
    """Сравнение букв: несколько разных по типу букв сразу, с теми же слоями и настройками.
    При смене настроек все они пересчитываются (по одной за такт, окно не замирает)."""

    CELL_W, CELL_H = 205, 245

    def __init__(self, main: "MainWindow"):
        super().__init__()
        self.main = main
        self._gen = 0
        self._queue: List[int] = []
        self.cells: List[tuple] = []  # (оболочка, заголовок, холст, символ, имя глифа)
        self._chars = ""
        v = QVBoxLayout(self)
        top = QHBoxLayout()
        top.addWidget(QLabel("Контрольные буквы:"))
        self.edit = QLineEdit("bdnagkxHRSeЖ")
        self.edit.setToolTip("Буквы, которые показываются одновременно. Подберите разные по типу: "
                             "с кривыми и только с прямыми, со стыками арки и стойки, с диагоналями.")
        self.edit.editingFinished.connect(self.refresh)
        btn = QPushButton("Обновить")
        btn.clicked.connect(lambda _=False: self.refresh())
        top.addWidget(self.edit, 1)
        top.addWidget(btn)
        v.addLayout(top)
        self.status = QLabel("")
        v.addWidget(self.status)
        self.scroll = QScrollArea()
        self.scroll.setWidgetResizable(True)
        self.body = QWidget()
        self.grid = QGridLayout(self.body)
        self.grid.setAlignment(Qt.AlignmentFlag.AlignTop | Qt.AlignmentFlag.AlignLeft)
        self.scroll.setWidget(self.body)
        v.addWidget(self.scroll, 1)

    def canvases(self) -> List[GlyphCanvas]:
        return [c[2] for c in self.cells]

    def _rebuild(self):
        store = self.main.store
        for wrap, _t, _c, _ch, _n in self.cells:
            wrap.setParent(None)
            wrap.deleteLater()
        self.cells = []
        seen = set()
        for ch in self.edit.text():
            name = store.name_of_char.get(ch) if store else None
            if not name or ch in seen:
                continue
            seen.add(ch)
            wrap = QWidget()
            lv = QVBoxLayout(wrap)
            lv.setContentsMargins(2, 2, 2, 2)
            title = QLabel(ch)
            cv = GlyphCanvas(interactive=False)
            cv.setFixedSize(self.CELL_W, self.CELL_H)
            cv.layers = self.main.canvas.layers  # одни и те же слои, что на основном холсте
            cv.color_mode = self.main.canvas.color_mode
            lv.addWidget(title)
            lv.addWidget(cv)
            self.cells.append((wrap, title, cv, ch, name))
        self._chars = self.edit.text()
        self._relayout()

    def _relayout(self):
        cols = max(1, (self.scroll.viewport().width() - 10) // (self.CELL_W + 12))
        for i, cell in enumerate(self.cells):
            self.grid.addWidget(cell[0], i // cols, i % cols)

    def resizeEvent(self, ev):
        super().resizeEvent(ev)
        self._relayout()

    def refresh(self):
        main = self.main
        if main.store is None:
            self.status.setText("Сначала откройте шрифт.")
            return
        if self.edit.text() != self._chars or not self.cells:
            self._rebuild()
        self._gen += 1
        self._queue = list(range(len(self.cells)))
        self.status.setText(f"Пересчёт: {len(self.cells)} букв…")
        gen = self._gen
        QTimer.singleShot(0, lambda: self._step(gen))

    def _step(self, gen: int):
        if gen != self._gen:
            return  # настройки изменились раньше, чем закончился пересчёт: он уже не нужен
        main = self.main
        if not self._queue:
            self.status.setText(f"Букв: {len(self.cells)}. Настройки применены к каждой по её классу.")
            return
        _wrap, title, cv, ch, name = self.cells[self._queue.pop(0)]
        try:
            raw = main.store.raw(name, main.profiles["curved"].remove_overlaps)
            a = analyze(raw, main.store.upm, main.profiles)
            cv.set_analysis(a, main.store.ascent, main.store.descent, True)
            c = a.counts()
            cover = a.scores.get("cover")
            title.setText(f"{ch}   {'прямые' if a.glyph_class == 'straight' else 'кривые'}   "
                          f"штрихов {c['s_line'] + c['s_arc']}   "
                          + (f"{cover * 100:.0f}%" if cover is not None else ""))
        except Exception as ex:
            title.setText(f"{ch}: ошибка {ex}")
        QTimer.singleShot(0, lambda: self._step(gen))


class ExportDialog(QDialog):
    """Параметры сохранения: оси или тонкое начертание в новый OTF."""

    def __init__(self, parent, mode: str, thin_ratio: float):
        super().__init__(parent)
        self.setWindowTitle("Сохранение в OTF")
        self.setMinimumWidth(480)
        form = QFormLayout(self)
        self.cb_mode = QComboBox()
        self.cb_mode.addItems(["Оси", "Тонкое начертание"])
        self.cb_mode.setCurrentIndex(0 if mode == "axes" else 1)
        self.cb_thin = QComboBox()
        self.cb_thin.addItems(["Доля от исходной толщины", "Одинаковая толщина всех штрихов"])
        self.sp_ratio = QDoubleSpinBox()
        self.sp_ratio.setRange(0.05, 1.0)
        self.sp_ratio.setSingleStep(0.05)
        self.sp_ratio.setDecimals(2)
        self.sp_ratio.setValue(thin_ratio)
        self.sp_abs = QDoubleSpinBox()
        self.sp_abs.setRange(1.0, 400.0)
        self.sp_abs.setValue(40.0)
        self.sp_abs.setSuffix(" ед.")
        self.sp_axis = QDoubleSpinBox()
        self.sp_axis.setRange(1.0, 100.0)
        self.sp_axis.setValue(10.0)
        self.sp_axis.setSuffix(" ед.")
        self.cb_axis_style = QComboBox()
        self.cb_axis_style.addItems(["Центральные линии (без расширения)", "Узкие залитые ленты (как раньше)"])
        self.cb_axis_style.setCurrentIndex(1)  # сохраняем прежнее поведение по умолчанию
        self.cb_axis_style.setToolTip(
            "Центральные линии записываются в CFF CharString без построения лент. "
            "Прямые оси остаются одиночными; кривые проходят назад по той же геометрии, "
            "чтобы автоматическое замыкание CFF не создавало хорду. "
            "Режим предназначен для редактирования путей в редакторе шрифтов."
        )
        self.sp_fit = QDoubleSpinBox()
        self.sp_fit.setRange(0.1, 5.0)
        self.sp_fit.setSingleStep(0.1)
        self.sp_fit.setValue(0.8)
        self.sp_fit.setSuffix(" ед.")
        self.cb_pathops = QCheckBox("Объединять наложения (skia-pathops)")
        self.cb_pathops.setChecked(HAVE_PATHOPS)
        self.cb_pathops.setEnabled(HAVE_PATHOPS)
        if not HAVE_PATHOPS:
            self.cb_pathops.setToolTip("Нужен пакет skia-pathops:  pip install skia-pathops\n"
                                       "Без него ленты штрихов остаются отдельными перекрывающимися контурами.")
        self.ed_suffix = QLineEdit(" Axes" if mode == "axes" else " Thin")
        self.ed_suffix.setToolTip("Добавляется к названию шрифта, чтобы файл не конфликтовал с исходным.")
        self.cb_keep = QCheckBox("Глифы без найденных штрихов оставить как в исходном шрифте")
        self.cb_keep.setChecked(True)
        form.addRow("Что сохранять:", self.cb_mode)
        form.addRow("Толщина тонкого начертания:", self.cb_thin)
        form.addRow("   доля от исходной:", self.sp_ratio)
        form.addRow("   одинаковая толщина штрихов:", self.sp_abs)
        form.addRow("Представление осей:", self.cb_axis_style)
        form.addRow("Ширина узкой ленты:", self.sp_axis)
        form.addRow("Допуск подгонки кривых:", self.sp_fit)
        form.addRow(self.cb_pathops)
        form.addRow("Добавить к названию шрифта:", self.ed_suffix)
        form.addRow(self.cb_keep)
        info = QLabel("Параметры применяются ко всем глифам. Центральные линии записываются без расширения: "
                      "это пути для редактора, а не корректное отображаемое начертание. "
                      "Из-за отсутствия открытых контуров в CFF прямые оси записываются одиночными, "
                      "а кривые возвращаются по себе назад, чтобы не создавать закрывающие хорды. "
                      "Поддерживается OTF с контурами CFF; ширины и метрики сохраняются.")
        info.setWordWrap(True)
        form.addRow(info)
        btns = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        btns.accepted.connect(self.accept)
        btns.rejected.connect(self.reject)
        form.addRow(btns)
        self.cb_mode.currentIndexChanged.connect(self._sync)
        self.cb_thin.currentIndexChanged.connect(self._sync)
        self.cb_axis_style.currentIndexChanged.connect(self._sync)
        self._sync()

    def _sync(self, *_):
        axes = self.cb_mode.currentIndex() == 0
        self.cb_thin.setEnabled(not axes)
        self.sp_ratio.setEnabled(not axes and self.cb_thin.currentIndex() == 0)
        self.sp_abs.setEnabled(not axes and self.cb_thin.currentIndex() == 1)
        centerline = self.cb_axis_style.currentIndex() == 0
        self.cb_axis_style.setEnabled(axes)
        self.sp_axis.setEnabled(axes and not centerline)
        self.sp_fit.setEnabled(not (axes and centerline))
        self.cb_pathops.setEnabled(HAVE_PATHOPS and not (axes and centerline))
        if self.ed_suffix.text().strip() in ("", "Axes", "Thin"):
            self.ed_suffix.setText(" Axes" if axes else " Thin")

    def options(self) -> ExportOptions:
        return ExportOptions(
            mode="axes" if self.cb_mode.currentIndex() == 0 else "thin",
            thin_mode="ratio" if self.cb_thin.currentIndex() == 0 else "absolute",
            thin_ratio=self.sp_ratio.value(),
            thin_abs=self.sp_abs.value(),
            axis_width=self.sp_axis.value(),
            axis_style="centerline" if self.cb_axis_style.currentIndex() == 0 else "ribbon",
            fit_tol=self.sp_fit.value(),
            use_pathops=self.cb_pathops.isChecked(),
            name_suffix=self.ed_suffix.text(),
            keep_unprocessed=self.cb_keep.isChecked(),
        )


class MainWindow(QMainWindow):
    def __init__(self, path: Optional[str] = None):
        super().__init__()
        self.setWindowTitle("StrokeLab: общая ось, предсказание и сохранение")
        self.resize(1440, 860)
        self.setAcceptDrops(True)
        self.store: Optional[FontStore] = None
        self.current: Optional[str] = None
        self.summary = ""
        self.model = GlyphModel()
        # два набора настроек: для букв, где все рёбра прямые, и для букв с кривыми
        self.profiles: Dict[str, Params] = {"straight": Params(), "curved": Params()}
        self.edit_key = "curved"  # какой набор сейчас показан в полях
        self._bound: List[Tuple[str, QWidget]] = []
        self._loading = False
        self.preview_dirty = True
        self._build_ui()
        QShortcut(QKeySequence("Ctrl+O"), self, activated=self.open_dialog)
        if path:
            QTimer.singleShot(0, lambda: self.load_font(path))
        else:
            QTimer.singleShot(0, self.open_dialog)

    # --- построение интерфейса -------------------------------------------------
    def _spin(self, field_name, lo, hi, step, suffix="", decimals=1) -> QDoubleSpinBox:
        sp = QDoubleSpinBox()
        sp.setRange(lo, hi)
        sp.setSingleStep(step)
        sp.setDecimals(decimals)
        sp.setValue(float(getattr(self.profiles[self.edit_key], field_name)))
        if suffix:
            sp.setSuffix(suffix)
        sp.setKeyboardTracking(False)  # пересчёт после Enter или потери фокуса, а не на каждую цифру
        sp.valueChanged.connect(self.on_param_changed)
        self._bound.append((field_name, sp))
        return sp

    def _check(self, field_name, text, tip="") -> QCheckBox:
        cb = QCheckBox(text)
        cb.setChecked(bool(getattr(self.profiles[self.edit_key], field_name)))
        if tip:
            cb.setToolTip(tip)
        cb.toggled.connect(self.on_param_changed)
        self._bound.append((field_name, cb))
        return cb

    def _add_rows(self, layout, rows):
        """rows: (поле Params, подпись, подсказка, min, max, шаг, суффикс, знаков)."""
        for field_name, text, tip, lo, hi, step, suffix, dec in rows:
            lab = QLabel(text)
            lab.setWordWrap(True)
            sp = self._spin(field_name, lo, hi, step, suffix, dec)
            lab.setToolTip(tip)
            sp.setToolTip(tip)
            layout.addWidget(lab)
            layout.addWidget(sp)

    def _build_ui(self):
        root = QWidget()
        self.setCentralWidget(root)
        h = QHBoxLayout(root)

        # слева: шрифт и список глифов
        left = QWidget()
        left.setFixedWidth(290)
        lv = QVBoxLayout(left)
        lv.setContentsMargins(0, 0, 0, 0)
        btn = QPushButton("Открыть шрифт…  (Ctrl+O)")
        btn.clicked.connect(self.open_dialog)
        self.font_label = QLabel("Шрифт не загружен")
        self.font_label.setWordWrap(True)
        self.filter_edit = QLineEdit()
        self.filter_edit.setPlaceholderText("Нужные буквы, например: Hbonaekx")
        self.filter_edit.textChanged.connect(self.on_filter)
        self.cb_by_name = QCheckBox("искать по имени глифа (uni0416, a.sc …)")
        self.cb_by_name.toggled.connect(self.on_filter)
        self.view = QListView()
        self.view.setModel(self.model)
        self.view.setViewMode(QListView.ViewMode.IconMode)
        self.view.setResizeMode(QListView.ResizeMode.Adjust)
        self.view.setMovement(QListView.Movement.Static)
        self.view.setUniformItemSizes(True)
        self.view.setWrapping(True)
        self.view.setIconSize(QSize(GlyphModel.ICON, GlyphModel.ICON))
        self.view.setGridSize(QSize(88, 86))
        self.view.setSelectionMode(QListView.SelectionMode.SingleSelection)
        self.view.setEditTriggers(QListView.EditTrigger.NoEditTriggers)
        self.view.selectionModel().currentChanged.connect(self.on_current_changed)
        lv.addWidget(btn)
        lv.addWidget(self.font_label)
        lv.addWidget(self.filter_edit)
        lv.addWidget(self.cb_by_name)
        lv.addWidget(self.view, 1)
        h.addWidget(left)

        # центр: вкладки «Буква» (основной холст) и «Сравнение букв» (несколько букв сразу)
        self.canvas = GlyphCanvas()
        self.canvas.picked.connect(self.on_picked)
        self.preview = PreviewPanel(self)
        self.tabs = QTabWidget()
        self.tabs.addTab(self.canvas, "Буква")
        self.tabs.addTab(self.preview, "Сравнение букв")
        self.tabs.currentChanged.connect(self.on_tab_changed)
        h.addWidget(self.tabs, 1)

        # справа: прокручиваемая панель настроек и под ней окно сведений
        controls = QWidget()
        controls.setMinimumWidth(0)
        controls.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred)
        cv = QVBoxLayout(controls)

        g_save = QGroupBox("Экспорт осей и начертания")
        gsv = QVBoxLayout(g_save)
        b_axes = QPushButton("Оси в OTF…")
        b_axes.clicked.connect(lambda _=False: self.export_dialog("axes"))
        b_ufo = QPushButton("Открытые оси в UFO…")
        b_ufo.clicked.connect(lambda _=False: self.export_vector_dialog("ufo"))
        b_ds = QPushButton("Открытые оси: Designspace + UFO…")
        b_ds.clicked.connect(lambda _=False: self.export_vector_dialog("designspace"))
        b_glyphs = QPushButton("Открытые оси в Glyphs 3…")
        b_glyphs.clicked.connect(lambda _=False: self.export_vector_dialog("glyphs"))
        b_svg = QPushButton("Оси в SVG…")
        b_svg.clicked.connect(self.export_svg_dialog)
        b_thin = QPushButton("Тонкое начертание в OTF…")
        b_thin.clicked.connect(lambda _=False: self.export_dialog("thin"))
        note_save = QLabel(
            "UFO, Designspace и Glyphs содержат настоящие открытые центральные линии "
            "без заливки и дублирования. Designspace создаётся с одним мастером; "
            "для интерполяции понадобятся дополнительные мастера."
        )
        note_save.setWordWrap(True)
        export_buttons = (b_axes, b_ufo, b_ds, b_glyphs, b_svg, b_thin)
        for button in export_buttons:
            button.setMinimumWidth(0)
            button.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Fixed)
            button.setStyleSheet("QPushButton { text-align: left; padding-left: 8px; }")
        for widget in (*export_buttons, note_save):
            gsv.addWidget(widget)
        cv.addWidget(g_save)

        g_layers = QGroupBox("Слои")
        gl = QVBoxLayout(g_layers)
        self.layer_boxes: Dict[str, QCheckBox] = {}
        for key, text in (
            ("fill", "Заливка глифа"),
            ("strokes", "Пятна штрихов (полупрозрачные)"),
            ("axes", "Оси штрихов (пунктир: продления и мосты)"),
            ("joints", "Соединения осей (T/L/Y/X/G1), концы"),
            ("pred_thin", "Предсказание: тонкое начертание"),
            ("clip_thin", "Обрезка по исходному контуру"),
            ("pred_orig", "Проверка: ленты при исходной толщине"),
            ("edges", "Рёбра контура"),
            ("leftover", "Остатки: части рёбер вне штрихов"),
            ("normals", "Нормали внутрь материала"),
            ("vertices", "Узлы и углы"),
            ("labels", "Номера рёбер и штрихов"),
            ("metrics", "Базовая линия, метрики"),
        ):
            cb = QCheckBox(text)
            cb.setChecked(self.canvas.layers[key])
            cb.toggled.connect(self.sync_layers)
            gl.addWidget(cb)
            self.layer_boxes[key] = cb
        self.color_combo = QComboBox()
        self.color_combo.addItems(["Раскраска рёбер: по типу", "Раскраска рёбер: по номеру"])
        self.color_combo.currentIndexChanged.connect(self.sync_layers)
        gl.addWidget(self.color_combo)
        fit_btn = QPushButton("Вписать глиф в окно")
        fit_btn.setMinimumWidth(0)
        fit_btn.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Fixed)
        fit_btn.setStyleSheet("QPushButton { text-align: left; padding-left: 8px; }")
        fit_btn.clicked.connect(self.canvas.fit_view)
        gl.addWidget(fit_btn)
        cv.addWidget(g_layers)

        g_par = QGroupBox("Разбиение на рёбра (общие настройки)")
        gp = QVBoxLayout(g_par)
        self._add_rows(gp, [
            ("smooth_deg", "Узел гладкий, если излом не больше:", "", 0.5, 45.0, 0.5, " °", 1),
            ("flat_tol", "Кривая считается прямой, если прогиб не больше:",
             "Рёбро-кривая, отклоняющаяся от хорды не больше, считается прямой.\n"
             "От этого зависит класс буквы (только прямые или есть кривые) и пары «прямая + дуга».",
             0.0, 50.0, 0.1, " ед.", 2),
        ])
        gp.addWidget(self._check("split_extrema", "Резать кривые в экстремумах"))
        gp.addWidget(self._check("split_inflections", "Резать кривые в перегибах"))
        cb_ovl = self._check("remove_overlaps", "Удалять наложения контуров")
        if not HAVE_PATHOPS:
            cb_ovl.setEnabled(False)
            cb_ovl.setToolTip("Нужен пакет skia-pathops:  pip install skia-pathops")
        gp.addWidget(cb_ovl)
        cv.addWidget(g_par)

        g_prof = QGroupBox("Набор настроек")
        gq = QVBoxLayout(g_prof)
        self.profile_combo = QComboBox()
        self.profile_combo.addItems(["Редактирую: буквы с кривыми", "Редактирую: буквы только с прямыми"])
        self.profile_combo.currentIndexChanged.connect(self.on_profile_switch)
        self.class_label = QLabel("Класс текущей буквы: нет")
        self.class_label.setWordWrap(True)
        self.class_label.setMinimumWidth(0)
        self.class_label.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred)
        note = QLabel("Для каждой буквы автоматически берётся набор по её классу: если все рёбра "
                      "прямые, то набор «только с прямыми», иначе «с кривыми». "
                      "Поля ниже относятся к тому набору, который выбран здесь.")
        note.setWordWrap(True)
        gq.addWidget(self.profile_combo)
        gq.addWidget(self.class_label)
        gq.addWidget(note)
        cv.addWidget(g_prof)

        g_str = QGroupBox("Поиск штрихов (выбранный набор)")
        gs = QVBoxLayout(g_str)
        self._add_rows(gs, [
            ("max_dev", "Макс. разворот нормалей двух сторон штриха:",
             "Больше: находит штрихи с сильным контрастом и косыми концами, но больше ложных пар.",
             5.0, 80.0, 1.0, " °", 1),
            ("disc_slack", "Строгость проверки пустого круга (0.5 – 1.0):",
             "Круг с диаметром между сторонами штриха не должен задевать чужие рёбра.\n"
             "Меньше: штрихи заходят глубже в стыки. Больше: стыки аккуратнее, но короче штрихи.",
             0.5, 1.0, 0.05, "", 2),
            ("merge_deg", "Допуск угла при склейке осей:",
             "Куски, разорванные перемычкой (стойка H, перекладина T), склеиваются в один штрих,\n"
             "если их оси продолжают друг друга с таким допуском.", 1.0, 90.0, 1.0, " °", 1),
            ("line_merge_deg", "Допуск угла при склейке осей-прямых:",
             "Прямые куски склеиваются в один штрих, если их оси почти продолжают друг друга. "
             "Для дуг действует допуск выше. Мал: диагональ не слипается со стойкой (Y, K).",
             1.0, 45.0, 1.0, " °", 1),
            ("gap_factor", "Макс. разрыв при склейке:", "Разрыв между кусками оси, в толщинах штриха.",
             0.0, 12.0, 0.1, " толщ.", 1),
            ("blob_ratio", "Короче этого считать пятном:",
             "Если штрих короче, чем столько толщин, это пятно (точка, засечка, стык), а не штрих.\n0: не применять.",
             0.0, 4.0, 0.1, " толщ.", 1),
            ("dot_ratio", "Точка: отдельная деталь не крупнее:",
             "Отдельная деталь (точка над i и j, двоеточие), чей размер не больше столько толщин штриха,\n"
             "считается точкой, а не набором штрихов. Тонкое начертание получает ту же точку, "
             "уменьшенную вокруг центра. 0: не искать точки.", 0.0, 4.0, 0.1, " толщ.", 1),
        ])
        gs.addWidget(self._check(
            "allow_mixed", "Разрешить пары «прямая + дуга»",
            "По умолчанию штрих строится только между двумя прямыми или двумя дугами. Включение возвращает пары\n"
            "«прямая + дуга»: помогает там, где стойка слита с чашей (g), но даёт лишние оси в других буквах."))
        gs.addWidget(self._check(
            "smooth_joints", "Сглаживать стыки осей",
            "В узле между кусками оси «усы» (касательные) ложатся на одну прямую с узлом."))
        self._add_rows(gs, [
            ("snap_deg", "Привязка касательной в стыке к горизонтали/вертикали:",
             "Если касательная в стыке отличается от горизонтали или вертикали не больше\n"
             "чем на этот угол, она ставится точно по ним (узел в экстремуме).", 0.0, 45.0, 1.0, " °", 1),
        ])
        cv.addWidget(g_str)

        g_join = QGroupBox("Соединение осей и предсказание")
        gj = QVBoxLayout(g_join)
        gj.addWidget(self._check(
            "join_axes", "Соединять оси в общую ось",
            "Выключено: показываются только оси отдельных штрихов."))
        gj.addWidget(self._check(
            "arm_continue", "Дуга входит в стойку по продолжению кривой чаши",
            "Чаша b, d, p, q: от вершины чаши ось идёт одной дугой окружности и садится на ось стойки под\n"
            "заданным углом, как её рисует человек. Чаша остаётся овалом, стойка прямой. Если вершины нет на этом\n"
            "штрихе или дуга сильно уходит от измеренной оси, берётся следующий способ (по середине контура)."))
        self._add_rows(gj, [
            ("arm_angle", "Угол входа дуги чаши в стойку:",
             "Чаша b, d, p, q: от вершины чаши ось идёт одной дугой окружности и садится на ось стойки\n"
             "под этим углом. Меньше: положе, и дуга садится ниже по стойке; больше: круче и выше.",
             10.0, 80.0, 1.0, "°", 0),
        ])
        gj.addWidget(self._check(
            "medial_joints", "Стыки по настоящей середине контура",
            "Запасной способ для дуги, входящей в стойку, и мост между дугами одной чаши у узла контура:\n"
            "ось ведётся по линии, равноудалённой от контура с двух сторон (учитывая вырез между стойкой\n"
            "и чашей). Выключено: прежние шаблонные кривые (дуга по касательной, середина прохода)."))
        self._add_rows(gj, [
            ("ext_factor", "Допуск на продление оси к стыку:",
             "Сколько непокрытого участка стыка разрешено продлить, в толщинах продлеваемого штриха\n"
             "(плюс полтолщины цели). Больше: достаёт дальше, но риск ложных соединений.",
             0.0, 4.0, 0.1, " толщ.", 1),
            ("min_support", "Опора в контуре: вогнутых углов у стыка не меньше:",
             "Настоящий стык (T, L, X) оставляет в контуре вогнутые углы. 0: соединять без такой опоры.",
             0.0, 6.0, 1.0, "", 0),
            ("term_factor", "Торец: перед концом контур ближе этой доли толщины:",
             "Если перед концом оси контур почти сразу (ближе этой доли толщины), конец настоящий (торец), "
             "и соединять его не с чем.", 0.0, 1.0, 0.05, " толщ.", 2),
            ("cap_factor", "Достраивать свободный конец до видимого торца, если он ближе:",
             "Свободный конец, перед которым виден торец, достраивается до него. 0: не достраивать.",
             0.0, 4.0, 0.1, " толщ.", 1),
            ("thin_ratio", "Толщина тонкого начертания (доля от исходной, общая):",
             "Предсказание тонкого начертания: толщина штриха умножается на это число.",
             0.05, 1.0, 0.05, "", 2),
        ])
        cv.addWidget(g_join)

        for group_box in (g_save, g_layers, g_par, g_prof, g_str, g_join):
            group_box.setMinimumWidth(0)
            group_box.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred)

        legend = QLabel(
            "<b>Рёбра</b><br>"
            '<span style="color:#2878dc">■</span> прямое &nbsp; '
            '<span style="color:#eb7814">■</span> выпуклая дуга &nbsp; '
            '<span style="color:#1ea064">■</span> вогнутая дуга<br>'
            '<span style="color:#e12828">●</span> выпуклый угол &nbsp; '
            '<span style="color:#14141a">■</span> вогнутый угол &nbsp; '
            '<span style="color:#9640c8">◆</span> перегиб<br><br>'
            "<b>Штрихи</b><br>"
            "Цветное полупрозрачное пятно: штрих между двумя сторонами. "
            "Тёмная линия с точками: его ось. У одного штриха один цвет.<br>"
            '<span style="color:#8c8c96">■</span> серое: пятно без оси (точка, засечка, стык)<br>'
            '<span style="color:#cd0078">┈┈</span> розовый пунктир: часть ребра, не вошедшая в штрих<br><br>'
            "Штрих строится только между двумя прямыми или двумя дугами; пара «прямая + дуга» "
            "считается ошибочной и осей не даёт.<br><br>"
            "<b>Соединения осей</b><br>"
            "▽ T: конец штриха упирается в ось другого &nbsp; □ L: угол двух штрихов<br>"
            "△ Y: три конца в одном узле &nbsp; ⊗ X: оси пересекаются &nbsp; ○ G1: плавное продолжение<br>"
            "Цифра у узла (в подписях): сколько вогнутых углов контура рядом, то есть опора в контуре.<br>"
            '<span style="color:#e08000">┈┈</span> оранжевый пунктир: продление оси, мост или достройка до торца<br>'
            '<span style="color:#eb7800">○</span> оранжевое кольцо: свободный конец (не торец и не соединён)<br><br>'
            "<b>Предсказание</b><br>"
            "Синее: тонкое начертание (ленты штрихов по осям). Розовый пунктир: те же ленты при "
            "исходной толщине; они должны ложиться на контур буквы, процент совпадения указан в сведениях.<br><br>"
            "Клик по пятну или по ребру показывает подробности."
        )
        legend.setWordWrap(True)
        cv.addWidget(legend)
        cv.addStretch(1)

        scroll = QScrollArea()
        scroll.setMinimumWidth(0)
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        scroll.setWidget(controls)

        self.info = QPlainTextEdit()
        self.info.setReadOnly(True)

        right = QWidget()
        right.setFixedWidth(350)
        rv = QVBoxLayout(right)
        rv.setContentsMargins(0, 0, 0, 0)
        rv.addWidget(scroll, 3)
        rv.addWidget(self.info, 2)
        h.addWidget(right)

    # --- шрифт и глифы ---------------------------------------------------------
    def open_dialog(self):
        path, _ = QFileDialog.getOpenFileName(self, "Открыть шрифт", "", "Шрифты (*.otf *.ttf);;Все файлы (*)")
        if path:
            self.load_font(path)

    def load_font(self, path: str):
        try:
            store = FontStore(path)
        except Exception as ex:
            QMessageBox.critical(self, "Не удалось открыть шрифт", f"{path}\n\n{ex}")
            return
        self.store = store
        self.filter_edit.blockSignals(True)
        self.filter_edit.clear()
        self.filter_edit.blockSignals(False)
        self.model.set_store(store)
        self.font_label.setText(f"{store.title}\nглифов: {len(store.display_order)}, UPM: {store.upm}")
        self.setWindowTitle(f"StrokeLab: {store.title}")
        self.current = None
        self._select_row(self._preferred_row())
        self.preview_dirty = True
        self.refresh_preview_if_visible()

    def _preferred_row(self) -> int:
        for ch in "HhOoНн":
            n = self.store.name_of_char.get(ch) if self.store else None
            if n in self.model.names:
                return self.model.names.index(n)
        return 0

    def _select_row(self, row: int):
        if self.model.rowCount() == 0:
            self.current = None
            self.canvas.set_analysis(None, 0, 0, False)
            self.info.setPlainText("Нет глифов по такому фильтру.")
            return
        idx = self.model.index(max(0, min(row, self.model.rowCount() - 1)), 0)
        self.view.setCurrentIndex(idx)
        self.view.scrollTo(idx)

    def on_filter(self, *_):
        self.model.set_filter(self.filter_edit.text(), self.cb_by_name.isChecked())
        self._select_row(0)

    def on_current_changed(self, cur, _prev):
        if not cur.isValid() or not self.store:
            return
        self.current = self.model.names[cur.row()]
        self.show_glyph(refit=True)

    SHARED_FIELDS = ("smooth_deg", "flat_tol", "split_extrema", "split_inflections", "remove_overlaps", "thin_ratio")

    def _widgets_to_profile(self):
        """Значения полей записываются в редактируемый набор; общие параметры копируются в оба."""
        prm = self.profiles[self.edit_key]
        for field_name, w in self._bound:
            v = w.value() if isinstance(w, QDoubleSpinBox) else w.isChecked()
            cur = getattr(prm, field_name)
            if isinstance(cur, bool):
                setattr(prm, field_name, bool(v))
            elif isinstance(cur, int):
                setattr(prm, field_name, int(round(v)))
            else:
                setattr(prm, field_name, float(v))
        for other in self.profiles.values():
            if other is not prm:
                for f in self.SHARED_FIELDS:
                    setattr(other, f, getattr(prm, f))

    def _profile_to_widgets(self):
        prm = self.profiles[self.edit_key]
        self._loading = True
        for field_name, w in self._bound:
            w.blockSignals(True)
            v = getattr(prm, field_name)
            if isinstance(w, QDoubleSpinBox):
                w.setValue(float(v))
            else:
                w.setChecked(bool(v))
            w.blockSignals(False)
        self._loading = False

    def on_param_changed(self, *_):
        if self._loading:
            return
        self._widgets_to_profile()
        self.reanalyze()
        self.preview_dirty = True
        self.refresh_preview_if_visible()

    def refresh_preview_if_visible(self):
        if self.store is not None and self.tabs.currentWidget() is self.preview and self.preview_dirty:
            self.preview_dirty = False
            self.preview.refresh()

    def on_tab_changed(self, _idx):
        self.refresh_preview_if_visible()

    def export_dialog(self, mode: str):
        if self.store is None:
            QMessageBox.information(self, "Сохранение", "Сначала откройте шрифт.")
            return
        dlg = ExportDialog(self, mode, self.profiles["curved"].thin_ratio)
        if dlg.exec() != QDialog.DialogCode.Accepted:
            return
        opts = dlg.options()
        base = os.path.splitext(os.path.basename(self.store.path))[0]
        default = os.path.join(os.path.dirname(self.store.path),
                               f"{base}-{'Axes' if opts.mode == 'axes' else 'Thin'}.otf")
        path, _ = QFileDialog.getSaveFileName(self, "Сохранить OTF", default, "OpenType (*.otf)")
        if not path:
            return
        n = len(self.store.display_order)
        prog = QProgressDialog("Обработка глифов…", "Отмена", 0, n, self)
        prog.setWindowModality(Qt.WindowModality.WindowModal)
        prog.setMinimumDuration(0)

        def cb(i: int, total: int, name: str) -> bool:
            prog.setValue(i)
            prog.setLabelText(f"Глиф {i + 1} из {total}: {name}")
            QApplication.processEvents()
            return not prog.wasCanceled()

        try:
            rep = export_font(self.store.path, path, self.profiles, opts, cb)
        except Exception as ex:
            prog.close()
            QMessageBox.critical(self, "Не удалось сохранить", f"{type(ex).__name__}: {ex}")
            return
        prog.setValue(n)
        if rep.cancelled:
            QMessageBox.information(self, "Сохранение", "Отменено, файл не записан.")
            return
        lines = [f"Сохранено: {path}", f"Заменено глифов: {rep.done} из {rep.total}"]
        if rep.skipped:
            lines.append(f"Без найденных штрихов: {len(rep.skipped)} ({', '.join(rep.skipped[:8])}"
                         + ("…)" if len(rep.skipped) > 8 else ")"))
        if rep.failed:
            lines.append(f"Ошибок: {len(rep.failed)}; первые: " + "; ".join(rep.failed[:3]))
        if rep.warnings:
            lines.append("Предупреждения: " + "; ".join(rep.warnings[:3]))
        if opts.mode == "axes" and opts.axis_style == "centerline":
            lines.append("Записаны центральные пути без расширения. Прямые оси остаются одиночными; криволинейные пути проходят назад по той же геометрии, потому что CFF/OTF автоматически замыкает под-контуры и иначе добавил бы хорду.")
        elif opts.use_pathops is False or not HAVE_PATHOPS:
            lines.append("Наложения лент не объединены (нужен skia-pathops): уберите их в редакторе шрифтов.")
        QMessageBox.information(self, "Сохранение", "\n".join(lines))

    def export_vector_dialog(self, kind: str):
        if self.store is None:
            QMessageBox.information(self, "Сохранение", "Сначала откройте шрифт.")
            return
        base = os.path.splitext(os.path.basename(self.store.path))[0]
        suffixes = {
            "ufo": ("UFO", "ufo", "UFO 3 (*.ufo);;Все файлы (*)", export_ufo),
            "designspace": ("Designspace", "designspace", "Designspace (*.designspace);;Все файлы (*)", export_designspace),
            "glyphs": ("Glyphs 3", "glyphs", "Glyphs 3 (*.glyphs);;Все файлы (*)", export_glyphs),
        }
        if kind not in suffixes:
            QMessageBox.critical(self, "Экспорт", f"Неизвестный формат: {kind}")
            return
        title, ext, file_filter, exporter = suffixes[kind]
        default = os.path.join(os.path.dirname(self.store.path), f"{base}-Centerlines.{ext}")
        path, _ = QFileDialog.getSaveFileName(self, f"Сохранить открытые оси в {title}", default, file_filter)
        if not path:
            return
        if not path.lower().endswith("." + ext):
            path += "." + ext
        count = len(self.store.display_order)
        prog = QProgressDialog(f"Экспорт осей в {title}…", "Отмена", 0, count, self)
        prog.setWindowModality(Qt.WindowModality.WindowModal)
        prog.setMinimumDuration(0)

        def cb(i: int, total: int, name: str) -> bool:
            prog.setValue(i)
            prog.setLabelText(f"Глиф {i + 1} из {total}: {name}")
            QApplication.processEvents()
            return not prog.wasCanceled()

        try:
            rep = exporter(self.store.path, path, self.profiles, progress=cb)
        except Exception as ex:
            prog.close()
            QMessageBox.critical(self, f"Не удалось сохранить {title}", f"{type(ex).__name__}: {ex}")
            return
        prog.setValue(count)
        if rep.cancelled:
            QMessageBox.information(self, "Экспорт", "Операция отменена, итоговый файл не создан.")
            return
        lines = [f"Сохранено: {path}", f"Глифов с осями: {rep.done} из {rep.total}"]
        if kind == "ufo":
            lines.append("Создан пакет UFO 3 с открытыми путями в GLIF.")
        elif kind == "designspace":
            companion = os.path.splitext(path)[0] + ".ufo"
            lines.append(f"UFO-мастер: {companion}")
            lines.append("Это один источник без интерполяционных осей; дополнительные мастера можно добавить позднее.")
        elif kind == "glyphs":
            lines.append("Создан текстовый файл Glyphs 3 с открытыми путями.")
        if rep.skipped:
            lines.append(f"Без найденных осей: {len(rep.skipped)}")
        if rep.failed:
            lines.append(f"Ошибок: {len(rep.failed)}; первые: " + "; ".join(rep.failed[:3]))
        if rep.warnings:
            lines.append("Примечания: " + "; ".join(rep.warnings[:3]))
        QMessageBox.information(self, f"Экспорт {title}", "\n".join(lines))

    def export_svg_dialog(self):
        if self.store is None:
            QMessageBox.information(self, "Сохранение", "Сначала откройте шрифт.")
            return
        base = os.path.splitext(os.path.basename(self.store.path))[0]
        default = os.path.join(os.path.dirname(self.store.path), f"{base}-Centerlines.svg")
        path, _ = QFileDialog.getSaveFileName(self, "Сохранить центральные линии в SVG", default,
                                              "SVG (*.svg);;Все файлы (*)")
        if not path:
            return
        if not path.lower().endswith(".svg"):
            path += ".svg"
        n = len(self.store.display_order)
        prog = QProgressDialog("Экспорт осей в SVG…", "Отмена", 0, n, self)
        prog.setWindowModality(Qt.WindowModality.WindowModal)
        prog.setMinimumDuration(0)

        def cb(i: int, total: int, name: str) -> bool:
            prog.setValue(i)
            prog.setLabelText(f"Глиф {i + 1} из {total}: {name}")
            QApplication.processEvents()
            return not prog.wasCanceled()

        try:
            rep = export_svg(self.store.path, path, self.profiles,
                             axis_width=max(2.0, self.store.upm / 80.0), progress=cb)
        except Exception as ex:
            prog.close()
            QMessageBox.critical(self, "Не удалось сохранить SVG", f"{type(ex).__name__}: {ex}")
            return
        prog.setValue(n)
        if rep.cancelled:
            QMessageBox.information(self, "Сохранение", "Отменено, SVG не записан.")
            return
        lines = [f"Сохранено: {path}", f"Глифов с осями: {rep.done} из {rep.total}"]
        if rep.skipped:
            lines.append(f"Без найденных осей: {len(rep.skipped)}")
        if rep.failed:
            lines.append(f"Ошибок: {len(rep.failed)}; первые: " + "; ".join(rep.failed[:3]))
        QMessageBox.information(self, "Экспорт SVG", "\n".join(lines))

    def on_profile_switch(self, idx: int):
        self._widgets_to_profile()  # то, что правили в прежнем наборе, сохраняется
        self.edit_key = "straight" if idx == 1 else "curved"
        self._profile_to_widgets()

    def reanalyze(self, *_):
        if self.store and self.current:
            self.show_glyph(refit=False)

    def show_glyph(self, refit: bool):
        try:
            raw = self.store.raw(self.current, self.profiles["curved"].remove_overlaps)
            a = analyze(raw, self.store.upm, self.profiles)
        except Exception as ex:
            self.info.setPlainText(f"Ошибка анализа глифа {self.current}:\n{ex}")
            self.canvas.set_analysis(None, 0, 0, False)
            return
        self.canvas.set_analysis(a, self.store.ascent, self.store.descent, refit)
        c = a.counts()
        left = sum(1 for e in a.edges if e.uncovered())
        cls_name = "только прямые рёбра" if a.glyph_class == "straight" else "есть кривые"
        self.class_label.setText(f"Класс текущей буквы: {cls_name}. Применён набор «"
                                 f"{'только с прямыми' if a.glyph_class == 'straight' else 'с кривыми'}».")
        self.summary = (
            f"Глиф: {self.store.describe(self.current)}\n"
            f"Класс: {cls_name}\n"
            f"Ширина: {a.advance:g}\n"
            f"Контуров: {c['contours']},  рёбер: {c['edges']} (прямых {c['lines']}, дуг {c['curves']})\n"
            f"Углы: выпуклых {c['convex']}, вогнутых {c['concave']}, перегибов {c['inflection']}\n"
            f"Штрихов: {c['s_line'] + c['s_arc']} (прямых {c['s_line']}, дуг {c['s_arc']}), "
            f"точек: {c['s_dot']}, пятен: {c['s_blob']}\n"
            f"Рёбер с участками вне штрихов: {left}\n"
            f"Соединений: T {c['j_T']}, L {c['j_L']}, Y {c['j_Y']}, X {c['j_X']}, G1 {c['j_G1']};  "
            f"свободных концов: {c['open_ends']}"
        )
        if a.scores:
            self.summary += (f"\nПроверка предсказания по контуру: совпало {a.scores['cover'] * 100:.0f}% контура, "
                             f"лишнего {a.scores['over'] * 100:.0f}%")
        if a.stroke_error:
            self.summary += f"\nОШИБКА поиска штрихов: {a.stroke_error}"
        if a.skeleton_error:
            self.summary += f"\nОШИБКА соединения осей: {a.skeleton_error}"
        self.on_picked(None, None)

    # --- слои и информация -------------------------------------------------------
    def sync_layers(self, *_):
        for k, cb in self.layer_boxes.items():
            self.canvas.layers[k] = cb.isChecked()
        self.canvas.color_mode = "index" if self.color_combo.currentIndex() == 1 else "type"
        self.canvas.update()
        for cv_ in self.preview.canvases():
            cv_.color_mode = self.canvas.color_mode
            cv_.update()

    def on_picked(self, e: Optional[Edge], s: Optional[Stroke]):
        text = self.summary
        if s is not None:
            a0, a1 = (s.axis[0][0], s.axis[-1][3]) if s.axis else (None, None)
            lines = ["", f"Штрих S{s.sid}: {STROKE_NAMES[s.kind]}"]
            if s.kind == "blob":
                lines.append(f"Толщина (по лучам): от {s.w_min:.1f} до {s.w_max:.1f}")
            else:
                lines.append(f"Толщина: мин {s.w_min:.1f}, средняя {s.w_mean:.1f}, макс {s.w_max:.1f}")
                lines.append(f"Длина оси: {s.length:.1f}" + ("  (замкнутая)" if s.closed else ""))
                if a0 and a1:
                    lines.append(f"Ось: ({a0[0]:.1f}, {a0[1]:.1f}) → ({a1[0]:.1f}, {a1[1]:.1f})")
            lines.append(f"Кусков: {len(s.pieces)};  рёбра: {', '.join(str(i) for i in s.edge_ids())}")
            text += "\n" + "\n".join(lines)
        if e is not None:
            if e.kind == "line":
                kind = "прямая"
            elif e.bulge == "convex":
                kind = "выпуклая дуга (центр кривизны в материале)"
            else:
                kind = "вогнутая дуга (центр кривизны снаружи)"
            a0 = math.degrees(math.atan2(e.t_in[1], e.t_in[0]))
            a1 = math.degrees(math.atan2(e.t_out[1], e.t_out[0]))
            n = e.normal(0.5)
            unc = e.uncovered()
            lines = [
                "",
                f"Ребро №{e.gid}  (контур {e.contour}, ребро {e.index})",
                f"Тип: {kind}",
                f"Длина: {e.length:.1f}",
                f"Направление касательной: {a0:.1f}° → {a1:.1f}°",
                f"Поворот: {e.turn:+.1f}°",
                f"Нормаль внутрь (в середине): ({n[0]:+.3f}, {n[1]:+.3f})",
                f"Начало: ({e.pts[0][0]:.1f}, {e.pts[0][1]:.1f})",
                f"Конец:  ({e.pts[3][0]:.1f}, {e.pts[3][1]:.1f})",
                "Вне штрихов: " + (", ".join(f"{t0:.2f}–{t1:.2f}" for t0, t1 in unc) if unc else "нет"),
            ]
            if e.v_start:
                lines.append(f"Узел в начале: {VERTEX_NAMES[e.v_start.kind]}, излом {e.v_start.angle:+.1f}°")
            if e.v_end:
                lines.append(f"Узел в конце: {VERTEX_NAMES[e.v_end.kind]}, излом {e.v_end.angle:+.1f}°")
            text += "\n" + "\n".join(lines)
        if e is None and s is None and self.summary:
            text += "\n\nКликните по пятну штриха или по ребру, чтобы увидеть подробности."
        self.info.setPlainText(text)

    # --- перетаскивание файла в окно ---------------------------------------------
    def dragEnterEvent(self, ev):
        if ev.mimeData().hasUrls():
            ev.acceptProposedAction()

    def dropEvent(self, ev):
        for u in ev.mimeData().urls():
            f = u.toLocalFile()
            if f.lower().endswith((".otf", ".ttf")):
                self.load_font(f)
                break


def main():
    app = QApplication(sys.argv)
    win = MainWindow(sys.argv[1] if len(sys.argv) > 1 else None)
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
