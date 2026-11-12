import pytest
import sympy as s

from hypra import BackendTimeoutError, NativeConfig, expand_epsilon

e=s.Symbol('epsilon')

def test_bounded_native_success():
    result=expand_epsilon(1/(1-e),e,2,config=NativeConfig(timeout_seconds=10,use_cache=False))
    assert result.polynomial==1+e+e**2

def test_bounded_native_timeout():
    with pytest.raises(BackendTimeoutError,match='exceeded'):
        expand_epsilon(1/(1-e),e,2,config=NativeConfig(timeout_seconds=0.000001,use_cache=False))

@pytest.mark.parametrize('value',[0,-1,True,float('inf'),float('nan'),'10'])
def test_timeout_validation(value):
    with pytest.raises(ValueError):
        NativeConfig(timeout_seconds=value)
