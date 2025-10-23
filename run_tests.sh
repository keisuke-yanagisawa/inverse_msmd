#!/bin/bash
# テスト実行スクリプト

set -e

echo "========================================"
echo "  inverse_msmd テストスイート"
echo "========================================"
echo ""

# 引数の処理
TEST_TYPE=${1:-all}

case $TEST_TYPE in
    unit)
        echo "単体テストを実行中..."
        pytest tests/unit -v -m "unit"
        ;;
    integration)
        echo "統合テストを実行中..."
        pytest tests/integration -v -m "integration"
        ;;
    fast)
        echo "高速テスト（visualとslowを除く）を実行中..."
        pytest -v -m "not slow and not visual"
        ;;
    all)
        echo "全てのテストを実行中..."
        pytest tests/ -v
        ;;
    coverage)
        echo "カバレッジレポート付きで全テストを実行中..."
        pytest tests/ -v --cov=inverse_msmd --cov-report=html --cov-report=term
        echo ""
        echo "HTMLカバレッジレポートが htmlcov/index.html に生成されました"
        ;;
    *)
        echo "使用方法: $0 [unit|integration|fast|all|coverage]"
        echo ""
        echo "  unit        - 単体テストのみ実行"
        echo "  integration - 統合テストのみ実行"
        echo "  fast        - 高速テスト（visualとslowマーカーを除く）"
        echo "  all         - 全てのテストを実行（デフォルト）"
        echo "  coverage    - カバレッジレポート付きで実行"
        exit 1
        ;;
esac

echo ""
echo "========================================"
echo "  テスト完了"
echo "========================================"