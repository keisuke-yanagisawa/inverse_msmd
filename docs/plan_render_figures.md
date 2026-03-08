# 実装計画: `render_figures` 引数によるPyMOL 3D描画の自動生成

## 概要

メインワークフロー（単体実行・バッチ処理）に `render_figures` フラグを追加し、
処理完了後に PyMOL による3D構造図を自動生成する機能を組み込む。

既存の `pymol_visualization.py` の関数群をそのまま活用する。
全出力ファイル（タンパク質PDB・リガンドSDF・プローブPDB）は
**プローブ座標系に統一済み**のため、追加の座標変換は不要。

---

## 変更対象ファイル

| ファイル | 変更内容 |
|---------|---------|
| `inverse_msmd/substructure_replacement.py` | `integrated_substructure_replacement` に `render_figures` 引数を追加 |
| `inverse_msmd/batch_processing.py` | `run_batch_processing` に `render_figures` 引数を追加 |
| `inverse_msmd/cli.py` | 両CLIに `--render-figures` オプションを追加 |

`pymol_visualization.py` 自体の変更は不要。

---

## 1. `substructure_replacement.py` の変更

### 1-1. 関数シグネチャ

```python
def integrated_substructure_replacement(
    ligand_file: str,
    protein_file: str,
    from_file: str,
    to_file: str,
    output_dir: str,
    match_index: Optional[int] = None,
    profile_dir: Optional[str] = None,
    probe_id: Optional[str] = None,
    csv_output: Optional[str] = None,
    image_output: Optional[str] = None,
    deduplicate_by_smiles: bool = False,
    skip_steric_clash_check: bool = False,
    render_figures: bool = False,          # ← 追加
) -> List[Dict[str, Union[str, float, int]]]:
```

### 1-2. docstring に追記（`skip_steric_clash_check` の後に）

```
    render_figures : bool, default=False
        3D構造図（PyMOL）を自動生成するかどうか。
        Trueの場合、各パターンの複合体図・統合図を出力ディレクトリに保存します。
        視点はプローブ分子のPCA（compute_probe_view）で自動計算されます。
        PyMOLがインストールされていない場合は警告を表示してスキップします。
```

### 1-3. 描画ロジックの挿入位置

関数末尾の `print(f"\n完了: ...")` の **直前**（現在の1144行目付近）に以下を挿入：

```python
    # 3D構造図の生成（オプション）
    if render_figures and results:
        print(f"\n3D構造図を生成中...")
        try:
            from .pymol_visualization import (
                compute_probe_view, render_complex, render_combined,
                render_probe_with_maps, _find_profile_files
            )

            probe_pdb = f"{to_file}.pdb"
            view = compute_probe_view(probe_pdb)

            # プロファイルファイルの検索
            profile_files = {}
            if calculate_scores and profile_dir and probe_id:
                profile_files = _find_profile_files(profile_dir, probe_id)

            # プローブ+マップ図（1回のみ、プロファイルがある場合）
            if profile_files:
                panel_b = str(output_path / "probe_map.png")
                render_probe_with_maps(
                    probe_pdb=probe_pdb,
                    profile_files=profile_files,
                    output_png=panel_b,
                    view=view,
                )
                print(f"  ✓ プローブ+マップ: {panel_b}")

            # 各パターンの描画
            for result in results:
                pat_idx = result['pattern_index']

                # 複合体図（タンパク質+リガンド）
                panel_a = str(output_path / f"pattern_{pat_idx}_complex.png")
                render_complex(
                    protein_pdb=result['protein_file'],
                    ligand_sdf=result['ligand_file'],
                    output_png=panel_a,
                    view=view,
                )
                print(f"  ✓ パターン {pat_idx} 複合体: {panel_a}")

                # 統合図（タンパク質+リガンド+プローブ+マップ）
                if profile_files:
                    panel_c = str(output_path / f"pattern_{pat_idx}_combined.png")
                    render_combined(
                        protein_pdb=result['protein_file'],
                        ligand_sdf=result['ligand_file'],
                        probe_pdb=probe_pdb,
                        profile_files=profile_files,
                        output_png=panel_c,
                        view=view,
                    )
                    print(f"  ✓ パターン {pat_idx} 統合図: {panel_c}")

            print(f"  ✓ 3D描画完了")
        except ImportError as e:
            print(f"  ⚠ 警告: PyMOLが利用できないため3D描画をスキップします: {e}")
        except Exception as e:
            print(f"  ⚠ 警告: 3D描画中にエラーが発生しました: {e}")
```

### 1-4. 生成されるファイル（単体実行時）

```
output_dir/
├── pattern_0_ligand_replaced.sdf      (既存)
├── pattern_0_protein_aligned.pdb      (既存)
├── pattern_0_complex.png              ← NEW: タンパク質+リガンド
├── pattern_0_combined.png             ← NEW: 統合図（プロファイルがある場合）
├── probe_map.png                      ← NEW: プローブ+マップ（プロファイルがある場合）
├── results.csv                        (既存、csv_output指定時)
└── ...
```

---

## 2. `batch_processing.py` の変更

### 2-1. `run_batch_processing` のシグネチャ

```python
def run_batch_processing(
    batch_csv: str,
    ligand_file: str,
    protein_file: str,
    probe_base_dir: str,
    profile_base_dir: str,
    output_base_dir: str,
    parallel: bool = False,
    max_workers: int = 4,
    continue_on_error: bool = True,
    log_file: Optional[str] = None,
    skip_steric_clash_check: bool = False,
    render_figures: bool = False,          # ← 追加
) -> BatchResult:
```

### 2-2. 描画ロジックの挿入位置

`return batch_result`（現在の853行目）の **直前** に以下を挿入：

```python
    # 3D構造図の生成（オプション）
    if render_figures and batch_result.num_success > 0:
        logger.info("3D構造図を生成中...")
        try:
            from inverse_msmd.pymol_visualization import render_batch_results as _render_batch

            # 成功したジョブを to_probe でグループ化
            from collections import defaultdict
            probe_groups = defaultdict(list)
            for job_result in batch_result.job_results:
                if job_result.status == "success":
                    probe_groups[job_result.to_probe].append(job_result.job_id)

            for to_probe, job_id_list in probe_groups.items():
                probe_pdb = str(Path(probe_base_dir) / f"{to_probe}.pdb")
                if not Path(probe_pdb).exists():
                    logger.warning(f"プローブPDBが見つかりません: {probe_pdb}")
                    continue

                logger.info(
                    f"プローブ {to_probe} のジョブ ({len(job_id_list)}件) を描画中..."
                )
                _render_batch(
                    output_base_dir=output_base_dir,
                    probe_pdb=probe_pdb,
                    profile_dir=profile_base_dir,
                    probe_id=to_probe,
                    job_ids=job_id_list,
                )

            logger.info("3D描画完了")
        except ImportError as e:
            logger.warning(f"PyMOLが利用できないため3D描画をスキップします: {e}")
        except Exception as e:
            logger.warning(f"3D描画中にエラーが発生しました: {e}")
```

### 2-3. 設計ポイント

- **to_probe ごとにグループ化**: バッチCSVでは各ジョブが異なる `to_probe` を持ちうるため、同じ `to_probe` のジョブをまとめて `render_batch_results` を呼ぶ。
- **`render_batch_results` の `job_ids` 引数**: 対象ジョブを限定するために使用。
- **Panel B（プローブ+マップ）** は `to_probe` ごとに1回だけ生成（`render_batch_results` 内部の設計）。
- **視点の自動計算**: `render_batch_results` 内部で `compute_probe_view` が `view=None` 時に自動呼び出しされるため、明示的な呼び出し不要。

---

## 3. `cli.py` の変更

### 3-1. `run_single_main` （inverse-msmd-run）

引数追加（`optional` グループ内、`--skip-steric-clash-check` の後）：

```python
    optional.add_argument(
        '--render-figures',
        action='store_true',
        help='3D構造図（PyMOL）を自動生成する'
    )
```

関数呼び出しに追加：

```python
        results = integrated_substructure_replacement(
            ...
            skip_steric_clash_check=args.skip_steric_clash_check,
            render_figures=args.render_figures,     # ← 追加
        )
```

設定表示に追加（既存の立体障害チェック表示の後）：

```python
    if args.render_figures:
        print(f"3D描画             : 有効")
```

### 3-2. `run_batch_main` （inverse-msmd-batch）

引数追加（既存の `--visualize` 等の後）：

```python
    optional.add_argument(
        '--render-figures',
        action='store_true',
        help='3D構造図（PyMOL）を自動生成する（各ジョブの複合体・統合図）'
    )
```

関数呼び出しに追加：

```python
        result = run_batch_processing(
            ...
            skip_steric_clash_check=args.skip_steric_clash_check,
            render_figures=args.render_figures,     # ← 追加
        )
```

> **注意**: 既存の `--visualize` は2D構造式のグリッド画像（matplotlib/RDKit）を生成する機能。
> `--render-figures` は3D構造図（PyMOL）を生成する機能で、役割が異なる。両方同時指定も可能。

---

## 4. 座標系と視点の整理

```
座標系: すべてプローブ(to_mol)座標系に統一済み
         ┌─────────────────────────────────────────┐
         │  integrated_substructure_replacement      │
         │                                           │
         │  タンパク質座標 ──(逆変換)──→ プローブ座標系  │
         │  リガンド座標   ──(逆変換)──→ プローブ座標系  │
         │  プローブ座標   ──(そのまま)→ プローブ座標系  │
         └─────────────────────────────────────────┘
                         ↓
         ┌─────────────────────────────────────────┐
         │  pymol_visualization（追加変換なし）        │
         │                                           │
         │  cmd.load(protein_aligned.pdb)  ← そのまま │
         │  cmd.load(ligand_replaced.sdf)  ← そのまま │
         │  cmd.load(probe.pdb)            ← そのまま │
         │                                           │
         │  視点: compute_probe_view(probe.pdb)       │
         │  → プローブのPCA主成分で決定               │
         │  → PC1,PC2が面内、PC3が法線               │
         │  → tilt_deg=45°で斜め上から見下ろし        │
         └─────────────────────────────────────────┘
```

---

## 5. エラーハンドリング方針

- PyMOLは**オプション依存**。`ImportError` を捕捉して警告のみ出し、メイン処理の結果は正常に返す。
- 描画中の他のエラー（ファイル不在等）も `Exception` で捕捉し、描画をスキップするだけ。
- メインの置換処理・スコア計算が成功していれば、描画の失敗で全体を失敗にしない。

---

## 6. 実装手順チェックリスト

1. [ ] `substructure_replacement.py`: シグネチャに `render_figures: bool = False` 追加
2. [ ] `substructure_replacement.py`: docstring に `render_figures` の説明を追加
3. [ ] `substructure_replacement.py`: 関数末尾に描画ロジックを挿入
4. [ ] `batch_processing.py`: `run_batch_processing` シグネチャに `render_figures` 追加
5. [ ] `batch_processing.py`: docstring 更新
6. [ ] `batch_processing.py`: `return batch_result` 直前に描画ロジック挿入
7. [ ] `cli.py`: `run_single_main` に `--render-figures` 引数追加 + 関数呼び出しに追加
8. [ ] `cli.py`: `run_batch_main` に `--render-figures` 引数追加 + 関数呼び出しに追加
9. [ ] テスト: PyMOLなし環境で `render_figures=False`（デフォルト）が既存動作を壊さないことを確認
10. [ ] テスト: `render_figures=True` + PyMOLなし環境で警告が出て正常終了することを確認
