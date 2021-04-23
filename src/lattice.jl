# Honeycomb lattice
Hhex2 = Lattice("Hhex2", [Point(PID(1), (0.0, 0.0)), Point(PID(2), (0.5, √3/6))],
                 vectors = [SVector(1.0, 0.0), SVector(0.5, √3/2)])

Hhex12 = Lattice("Hhex12", [Point(PID(1), (0.0, 0.0)), Point(PID(2), (1/2, √3/6)),
                            Point(PID(3), (1.0, 0.0)), Point(PID(4), (1.5, √3/6)),
                            Point(PID(5), (2.0, 0.0)), Point(PID(6), (2.5, √3/6)),
                            Point(PID(7), (0.5, √3/2)), Point(PID(8), (1.0, 2 *√3/3)),
                            Point(PID(9), (1.5, √3/2)), Point(PID(10), (2.0, 2 *√3/3)),
                            Point(PID(11), (2.5, √3/2)), Point(PID(12), (3.0, 2 *√3/3))],
                 vectors = [SVector(3.0, 0.0), SVector(1.0, √3)])

Hhex18 = Lattice("Hhex18", [Point(PID(1), (0.0, 0.0)), Point(PID(2), (1/2, √3/6)),
                            Point(PID(3), (1.0, 0.0)), Point(PID(4), (1.5, √3/6)),
                            Point(PID(5), (2.0, 0.0)), Point(PID(6), (2.5, √3/6)),
                            Point(PID(7), (0.5, √3/2)), Point(PID(8), (1.0, 2 *√3/3)),
                            Point(PID(9), (1.5, √3/2)), Point(PID(10), (2.0, 2 *√3/3)),
                            Point(PID(11), (2.5, √3/2)), Point(PID(12), (3.0, 2 *√3/3)),
                            Point(PID(13), (1.0, √3)), Point(PID(14), (1.5, 7*√3/6 )),
                            Point(PID(15), (2.0, √3)), Point(PID(16), (2.5, 7*√3/6)),
                            Point(PID(17), (3.0, √3)), Point(PID(18), (3.5, 7*√3/6))],
                 vectors = [SVector(3.0, 0.0), SVector(1.5, 1.5 *√3)])

Hhex24 = Lattice("Hhex24", [Point(PID(1), (0.0, 0.0)), Point(PID(2), (1/2, √3/6)),
                            Point(PID(3), (1.0, 0.0)), Point(PID(4), (1.5, √3/6)),
                            Point(PID(5), (2.0, 0.0)), Point(PID(6), (2.5, √3/6)),
                            Point(PID(7), (0.5, √3/2)), Point(PID(8), (1.0, 2 *√3/3)),
                            Point(PID(9), (1.5, √3/2)), Point(PID(10), (2.0, 2 *√3/3)),
                            Point(PID(11), (2.5, √3/2)), Point(PID(12), (3.0, 2 *√3/3)),
                            Point(PID(13), (1.0, √3)), Point(PID(14), (1.5, 7*√3/6 )),
                            Point(PID(15), (2.0, √3)), Point(PID(16), (2.5, 7*√3/6)),
                            Point(PID(17), (3.0, √3)), Point(PID(18), (3.5, 7*√3/6)),
                            Point(PID(19), (1.5, 3*√3/2)), Point(PID(20), (2.0, √3/6 + 3*√3/2)),
                            Point(PID(21), (2.5, 3*√3/2)), Point(PID(22), (3.0, √3/6 + 3*√3/2)),
                            Point(PID(23), (3.5, 3*√3/2)), Point(PID(24), (4.0, √3/6 + 3*√3/2))],
                 vectors = [SVector(3.0, 0.0), SVector(2.0, 2.0 *√3)])

H24 = Lattice("Hhex24", [Point(PID(1), (0.0, 0.0)), Point(PID(2), (0.5, -√3/6)),
                            Point(PID(3), (1.0, 0.0)), Point(PID(4), (3/2, -√3/6)),
                            Point(PID(5), (2.0, 0.0)), Point(PID(6), (-1/2, √3/2)),
                            Point(PID(7), (0.0, √3/3)), Point(PID(8), (1/2, √3/2)),
                            Point(PID(9), (1.0, √3/3)), Point(PID(10), (1.5, √3/2)),
                            Point(PID(11), (2.0, √3/3)), Point(PID(12), (2.5, √3/2)),
                            Point(PID(13), (-1/2, √3/2 + √3/3)), Point(PID(14), (0.0, √3)),
                            Point(PID(15), (0.5, √3/2 + √3/3)), Point(PID(16), (1.0, √3)),
                            Point(PID(17), (1.5, √3/2 + √3/3)), Point(PID(18), (2.0, √3)),
                            Point(PID(19), (2.5, √3/2 + √3/3)), Point(PID(20), (0.0, √3 + √3/3)),
                            Point(PID(21), (0.5, √3 + √3/2)), Point(PID(22), (1.0, √3 + √3/3)),
                            Point(PID(23), (1.5, √3 + √3/2)), Point(PID(24), (2.0, √3 + √3/3))],
                 vectors = [SVector(3.0, √3), SVector(0.0, 2.0 *√3)])
